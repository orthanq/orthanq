use crate::calling::haplotypes::{
    AlleleFreqDist, CandidateMatrix, PriorTypes, VariantCalls, VariantStatus,
};
use bio::stats::probs::adaptive_integration;
use bio::stats::{bayesian::model, LogProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_new::new;
use ordered_float::NotNan;
use statrs::function::beta::ln_beta;
use std::collections::HashMap;
use std::mem;

pub(crate) type AlleleFreq = NotNan<f64>;

#[derive(Hash, PartialEq, Eq, Clone, Debug, Derefable)]
pub(crate) struct HaplotypeFractions(#[deref] Vec<AlleleFreq>);

#[derive(Debug, new)]
pub(crate) struct Marginal {
    n_haplotypes: usize,
    upper_bond: NotNan<f64>,
    prior_info: PriorTypes,
}

impl Marginal {
    pub(crate) fn calc_marginal<
        F: FnMut(&<Self as model::Marginal>::Event, &<Self as model::Marginal>::Data) -> LogProb,
    >(
        &self,
        data: &Data,
        haplotype_index: usize,
        fractions: &mut Vec<AlleleFreq>,
        joint_prob: &mut F,
    ) -> LogProb {
        if haplotype_index == self.n_haplotypes {
            let event = HaplotypeFractions(fractions.to_vec());
            joint_prob(&event, data)
        } else {
            let fraction_upper_bound = self.upper_bond - fractions.iter().sum::<NotNan<f64>>();
            let mut density = |fraction| {
                let mut fractions = fractions.clone();
                fractions.push(fraction);
                self.calc_marginal(data, haplotype_index + 1, &mut fractions, joint_prob)
            };
            if haplotype_index == self.n_haplotypes - 1 {
                let second_if = density(fraction_upper_bound);
                second_if
            } else {
                if fraction_upper_bound == NotNan::new(0.0).unwrap() {
                    let last_if = density(NotNan::new(0.0).unwrap());
                    last_if
                } else {
                    //check prior info
                    if self.prior_info == PriorTypes::Diploid {
                        //sum 0.0, 0.5 and 1.0
                        let mut probs = Vec::new();
                        let mut diploid_points = |point, probs: &mut Vec<_>| {
                            let mut fractions = fractions.clone();
                            if fractions.iter().sum::<NotNan<f64>>() + point
                                <= NotNan::new(1.0).unwrap()
                            {
                                // this check is necessary to avoid combinations that sum up to more than 1.0.
                                probs.push(density(point));
                            } else {
                                ()
                            }
                        };
                        diploid_points(NotNan::new(0.0).unwrap(), &mut probs);
                        diploid_points(NotNan::new(0.5).unwrap(), &mut probs);
                        diploid_points(NotNan::new(1.0).unwrap(), &mut probs);
                        LogProb::ln_sum_exp(&probs)
                    } else if self.prior_info == PriorTypes::Uniform
                        || self.prior_info == PriorTypes::DiploidSubclonal
                    {
                        adaptive_integration::ln_integrate_exp(
                            density,
                            NotNan::new(0.0).unwrap(),
                            fraction_upper_bound,
                            NotNan::new(0.1).unwrap(),
                        )
                    } else {
                        panic!("uniform, prior or diploid-subclonal must be selected")
                    }
                }
            }
        }
    }
}

impl model::Marginal for Marginal {
    type Event = HaplotypeFractions;
    type Data = Data;
    type BaseEvent = HaplotypeFractions;

    fn compute<F: FnMut(&Self::Event, &Self::Data) -> LogProb>(
        &self,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let mut fractions: Vec<AlleleFreq> = Vec::new();
        self.calc_marginal(data, 0, &mut fractions, joint_prob)
    }
}

#[derive(Debug, new)]
pub(crate) struct Data {
    pub candidate_matrix: CandidateMatrix,
    pub variant_calls: VariantCalls,
    // pub kallisto_estimates: Vec<KallistoEstimate>
}

#[derive(Debug, new)]
pub(crate) struct Likelihood;

impl model::Likelihood<Cache> for Likelihood {
    type Event = HaplotypeFractions;
    type Data = Data;

    fn compute(&self, event: &Self::Event, data: &Self::Data, payload: &mut Cache) -> LogProb {
        // self.compute_kallisto(event, data, payload) + //comment it out for now
        self.compute_varlociraptor(event, data, payload)
    }
}

impl Likelihood {
    fn compute_kallisto(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
        _cache: &mut Cache,
    ) -> LogProb {
        // TODO compute likelihood using neg_binom on the counts and dispersion
        // in the data and the fractions in the events.
        //Later: use the cache to avoid redundant computations.
        // event
        //     .iter()
        //     .zip(data.kallisto_estimates.iter())
        //     .map(|(fraction, estimate)| {
        //         dbg!(&estimate, &fraction);
        //         neg_binom(
        //             *estimate.count,
        //             NotNan::into_inner(*fraction),
        //             *estimate.dispersion,
        //         )
        //     })
        //     .sum();
        LogProb::ln_one()
    }

    fn compute_varlociraptor(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
        _cache: &mut Cache,
    ) -> LogProb {
        let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
            data.candidate_matrix.values().cloned().collect();
        let variant_calls: Vec<AlleleFreqDist> = data
            .variant_calls
            .iter()
            .map(|(_, (_, afd))| afd.clone())
            .collect();
        let mut final_prob = LogProb::ln_one();
        candidate_matrix_values
            .iter()
            .zip(variant_calls.iter())
            .for_each(|((genotypes, covered), afd)| {
                let mut denom = NotNan::new(1.0).unwrap();
                let mut vaf_sum = NotNan::new(0.0).unwrap();
                event.iter().enumerate().for_each(|(i, fraction)| {
                    if genotypes[i] == VariantStatus::Present && covered[i as u64] {
                        vaf_sum += *fraction;
                    } else if genotypes[i] == VariantStatus::Unknown {
                        ()
                    } else if covered[i as u64] {
                        ()
                    } else if genotypes[i] == VariantStatus::NotPresent
                        && covered[i as u64] == false
                    {
                        denom -= *fraction;
                    }
                });
                if denom > NotNan::new(0.0).unwrap() {
                    vaf_sum /= denom;
                }
                //to overcome a bug that results in larger than 1.0 VAF. After around 10 - 15th decimal place, the value becomes larger.
                //In any case, for a direct query to the AFD VAFs (they contain 2 decimal places).
                vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round()).unwrap()
                    / NotNan::new(100.0).unwrap();
                if !afd.is_empty() {
                    final_prob += afd.vaf_query(&vaf_sum).unwrap();
                } else {
                    final_prob += LogProb::ln_one();
                }
            });
        final_prob
        //LogProb::ln_one()
    }
}

#[derive(Debug, new)]
pub(crate) struct Prior {
    prior: PriorTypes,
}

impl model::Prior for Prior {
    type Event = HaplotypeFractions;

    fn compute(&self, event: &Self::Event) -> LogProb {
        if self.prior == PriorTypes::Diploid {
            let mut prior_prob = LogProb::ln_one();
            event.iter().for_each(|fraction| {
                if *fraction == NotNan::new(0.0).unwrap() {
                    prior_prob += LogProb::from(Prob(1.0 / 3.0))
                } else if *fraction == NotNan::new(0.5).unwrap() {
                    prior_prob += LogProb::from(Prob(1.0 / 3.0))
                } else if *fraction == NotNan::new(1.0).unwrap() {
                    prior_prob += LogProb::from(Prob(1.0 / 3.0))
                } else {
                    prior_prob += LogProb::ln_zero()
                }
            });
            prior_prob
        } else if self.prior == PriorTypes::DiploidSubclonal {
            //diploid subclonal prior: don't allow for more than 4 fractions bearing greater than 0.0
            if event
                .iter()
                .filter(|&n| n > &NotNan::new(0.0).unwrap())
                .count()
                > 4
            {
                LogProb::ln_zero()
            } else {
                LogProb::ln_one()
            }
        } else {
            LogProb::ln_one()
        }
    }
}

#[derive(Debug, new)]
pub(crate) struct Posterior;

impl model::Posterior for Posterior {
    type Event = HaplotypeFractions;

    type BaseEvent = HaplotypeFractions;

    type Data = Data;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        // joint_prob calculates the joint probability from likelihood and prior
        joint_prob(event, data)
    }
}

#[derive(Debug, Derefable, Default)]
pub(crate) struct Cache(#[deref] HashMap<usize, HashMap<AlleleFreq, LogProb>>);

// fn recursive_vaf_query(
//     haplotype_index: usize,
//     fractions: &Vec<AlleleFreq>,
//     upper_bond: &NotNan<f64>,
//     afd: &AlleleFreqDist,
// ) -> LogProb {
//     let n_haplotypes = 1;
//     let mut vaf_sum = NotNan::new(0.0).unwrap();
//     if haplotype_index == n_haplotypes {
//         vaf_sum += *fractions[0];
//         if !afd.is_empty() {
//             afd.vaf_query(&vaf_sum).unwrap()
//         } else {
//             LogProb::ln_one()
//         }
//     } else {
//         let fraction_upper_bound = *upper_bond - fractions.iter().sum::<NotNan<f64>>();
//         let mut density = |fraction| {
//             let mut fractions = fractions.clone();
//             fractions.push(fraction);
//             recursive_vaf_query(haplotype_index + 1, &mut fractions, upper_bond, afd)
//         };
//         if fraction_upper_bound == NotNan::new(0.0).unwrap() {
//             density(NotNan::new(0.0).unwrap())
//         } else {
//             adaptive_integration::ln_integrate_exp(
//                 density,
//                 NotNan::new(0.0).unwrap(),
//                 fraction_upper_bound,
//                 NotNan::new(0.1).unwrap(),
//             )
//         }
//     }
// }

pub(crate) fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);
    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);
    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    LogProb((p1 - b + p2) - (x + n).ln())
}
