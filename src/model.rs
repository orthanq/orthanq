use crate::calling::haplotypes::{AlleleFreqDist, CandidateMatrix, VariantCalls};
use bio::stats::probs::adaptive_integration;
use bio::stats::{bayesian::model, LogProb};
use bv::BitVec;
use derefable::Derefable;
use derive_new::new;
use ordered_float::NotNan;
use std::collections::HashMap;

pub(crate) type AlleleFreq = NotNan<f64>;

#[derive(Hash, PartialEq, Eq, Clone, Debug, Derefable)]
pub(crate) struct HaplotypeFractions(#[deref] Vec<AlleleFreq>);

#[derive(Debug, new)]
pub(crate) struct Marginal {
    n_haplotypes: usize,
    upper_bond: NotNan<f64>,
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
                density(fraction_upper_bound)
            } else {
                if fraction_upper_bound == NotNan::new(0.0).unwrap() {
                    density(NotNan::new(0.0).unwrap())
                } else {
                    adaptive_integration::ln_integrate_exp(
                        density,
                        NotNan::new(0.0).unwrap(),
                        fraction_upper_bound,
                        NotNan::new(0.1).unwrap(),
                    )
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
}

#[derive(Debug, new)]
pub(crate) struct Likelihood {
    normalization: bool,
}

impl model::Likelihood<Cache> for Likelihood {
    type Event = HaplotypeFractions;
    type Data = Data;

    fn compute(&self, event: &Self::Event, data: &Self::Data, payload: &mut Cache) -> LogProb {
        self.compute_varlociraptor(event, data, payload)
    }
}

impl Likelihood {
    fn compute_varlociraptor(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
        _cache: &mut Cache,
    ) -> LogProb {
        let candidate_matrix: Vec<(Vec, BitVec)> =
            data.candidate_matrix.values().cloned().collect();
        let variant_calls: Vec<AlleleFreqDist> = data
            .variant_calls
            .iter()
            .map(|(_, afd)| afd.clone())
            .collect();
        if self.normalization {
            candidate_matrix
                .iter()
                .zip(variant_calls.iter())
                .map(|((genotypes, covered), afd)| {
                    let mut denom = NotNan::new(1.0).unwrap();
                    let mut vaf_sum = NotNan::new(0.0).unwrap();
                    event.iter().enumerate().for_each(|(i, fraction)| {
                        if genotypes[i as u64] && covered[i as u64] {
                            vaf_sum += *fraction;
                        } else if covered[i as u64] {
                            ()
                        } else {
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
                        afd.vaf_query(&vaf_sum).unwrap()
                    } else {
                        LogProb::ln_one()
                    }
                })
                .sum()
        } else {
            candidate_matrix
                .iter()
                .zip(variant_calls.iter())
                .map(|((genotypes, covered), afd)| {
                    let mut vaf_sum = NotNan::new(0.0).unwrap();
                    event.iter().enumerate().for_each(|(i, fraction)| {
                        if genotypes[i as u64] && covered[i as u64] {
                            vaf_sum += *fraction;
                        } else if covered[i as u64] {
                            ()
                        }
                    });
                    if !afd.is_empty() {
                        afd.vaf_query(&vaf_sum).unwrap()
                    } else {
                        LogProb::ln_one()
                    }
                })
                .sum()
        }
        //LogProb::ln_one()
    }
}

#[derive(Debug, new)]
pub(crate) struct Prior;

impl model::Prior for Prior {
    type Event = HaplotypeFractions;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        // flat prior for now
        LogProb::ln_one()
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
