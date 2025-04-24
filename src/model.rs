use crate::calling::haplotypes::haplotypes::{
    AlleleFreqDist, CandidateMatrix, Haplotype, HaplotypeGraph, PriorTypes, VariantCalls,
};

use bio::stats::probs::adaptive_integration;
use bio::stats::{bayesian::model, LogProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_new::new;
use ordered_float::NotNan;
use petgraph::visit::Bfs;
use std::collections::HashMap;
pub type AlleleFreq = NotNan<f64>;

#[derive(Hash, PartialEq, Eq, Clone, Debug, Derefable, PartialOrd)]
pub struct HaplotypeFractions(#[deref] pub Vec<AlleleFreq>);

#[derive(Debug, new)]
pub(crate) struct Marginal {
    n_haplotypes: usize,
    haplotypes: Vec<Haplotype>,
    prior_info: PriorTypes,
    haplotype_graph: Option<HaplotypeGraph>,
    enable_equivalence_class_constraint: bool,
    application: String,
}

impl Marginal {
    pub(crate) fn calc_marginal<
        F: FnMut(&<Self as model::Marginal>::Event, &<Self as model::Marginal>::Data) -> LogProb,
    >(
        &self,
        data: &Data,
        haplotype_index: usize,
        fractions: &mut [AlleleFreq],
        joint_prob: &mut F,
    ) -> LogProb {
        if haplotype_index == self.n_haplotypes {
            let event = HaplotypeFractions(fractions.to_vec());
            joint_prob(&event, data)
        } else {
            let fraction_upper_bound =
                NotNan::new(1.00).unwrap() - fractions.iter().sum::<NotNan<f64>>();
            let mut density = |fraction| {
                if self.enable_equivalence_class_constraint
                    && fraction > NotNan::new(0.0).unwrap()
                    && fractions.len() > 1
                {
                    if self.application == "hla".to_string() {
                        // only if the fraction for the current has greater than 0.0
                        let current_haplotype = &self.haplotypes[haplotype_index];
                        let splitted = &self.haplotypes[haplotype_index]
                            .split(':')
                            .collect::<Vec<&str>>();
                        let haplotype_group =
                            Haplotype(splitted[0].to_owned() + &":" + splitted[1]);

                        //find the index of the (haplotype, haplotype_group) in graph
                        if let Some(haplotype_graph) = &self.haplotype_graph {
                            // query node index
                            let index = haplotype_graph
                                .get_node_index(&(current_haplotype.clone(), haplotype_group))
                                .unwrap();
                            // step through the graph and sum incoming edges into the node weight
                            let mut bfs = Bfs::new(&**haplotype_graph, index);

                            while let Some(nx) = bfs.next(&**haplotype_graph) {
                                // we can access `graph` mutably here still
                                let haplotype_query = &haplotype_graph[nx].0;
                                for (h, f) in self.haplotypes[0..haplotype_index]
                                    .to_vec()
                                    .iter()
                                    .zip(fractions[0..haplotype_index].to_vec().iter())
                                {
                                    if (h == haplotype_query) && (f > &NotNan::new(0.0).unwrap()) {
                                        return LogProb::ln_zero();
                                    }
                                }
                            }
                        }
                    }
                    // else if self.application == "virus".to_string() {
                    //TODO: explore other methods.
                    // }
                }
                // dbg!(&fractions);
                let mut fractions = fractions.to_vec();
                fractions.push(fraction);
                self.calc_marginal(data, haplotype_index + 1, &mut fractions, joint_prob)
            };

            if haplotype_index == self.n_haplotypes - 1 {
                density(fraction_upper_bound)
            } else {
                if fraction_upper_bound == NotNan::new(0.0).unwrap() {
                    density(NotNan::new(0.0).unwrap())
                } else {
                    //check prior info
                    if self.prior_info == PriorTypes::Diploid {
                        //sum 0.0, 0.5 and 1.0
                        let mut probs = Vec::new();
                        let mut diploid_points = |point, probs: &mut Vec<_>| {
                            let fractions = fractions.to_vec();
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
pub struct Data {
    pub candidate_matrix: CandidateMatrix,
    pub variant_calls: VariantCalls,
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
    fn compute_varlociraptor(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
        _cache: &mut Cache,
    ) -> LogProb {
        let candidate_matrix_values: Vec<(BitVec, BitVec)> =
            data.candidate_matrix.values().cloned().collect();
        let variant_calls: Vec<(AlleleFreqDist, i32)> = data
            .variant_calls
            .iter()
            .map(|(_, (_,_, _, afd, coverage))| (afd.clone(), *coverage))
            .collect();
        let mut final_prob = LogProb::ln_one();
        candidate_matrix_values
            .iter()
            .zip(variant_calls.iter())
            .for_each(|((genotypes, covered), (afd, cov))| {
                if *cov != 0 {
                    let mut denom = NotNan::new(1.0).unwrap();
                    let mut vaf_sum = NotNan::new(0.0).unwrap();
                    event.iter().enumerate().for_each(|(i, fraction)| {
                        if genotypes[i as u64] && covered[i as u64] {
                            vaf_sum += *fraction;
                        }
                        // else if genotypes[i] == VariantStatus::Unknown || covered[i as u64] {
                        //     ()
                        // }
                        else if !genotypes[i as u64] && !covered[i as u64] {
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
                if *fraction == NotNan::new(0.0).unwrap()
                    || *fraction == NotNan::new(0.5).unwrap()
                    || *fraction == NotNan::new(1.0).unwrap()
                {
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
