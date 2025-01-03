use crate::calling::haplotypes::haplotypes::SimilarL;
use crate::calling::haplotypes::haplotypes::{
    AlleleFreqDist, CandidateMatrix, Haplotype, HaplotypeGraph, PriorTypes, VariantCalls,
    VariantStatus,
};

use bio::stats::probs::adaptive_integration;
use bio::stats::{bayesian::model, LogProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_new::new;
use ordered_float::NotNan;
use petgraph::visit::Bfs;
use petgraph::Graph;
use petgraph::Undirected;
use std::collections::{BTreeMap, BTreeSet, HashMap};
pub type AlleleFreq = NotNan<f64>;

#[derive(Hash, PartialEq, Eq, Clone, Debug, Derefable, PartialOrd)]
pub struct HaplotypeFractions(#[deref] pub Vec<AlleleFreq>);

#[derive(Debug, new)]
pub(crate) struct Marginal {
    n_haplotypes: usize,
    haplotypes: Vec<Haplotype>,
    upper_bond: NotNan<f64>,
    prior_info: PriorTypes,
    haplotype_graph: Option<HaplotypeGraph>,
    // similar_haplotypes: Option<BTreeMap<Haplotype, BTreeSet<Haplotype>>>,
    // similarl_map: Option<SimilarL>,
    distance_threshold: Option<usize>,
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
            let fraction_upper_bound = self.upper_bond - fractions.iter().sum::<NotNan<f64>>();
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
                        let mut haplotype_group =
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
                        
                        // updated approach: 
                        // // new approach:Let L be the set of haplotypes predicted by the LP.
                        // // Let R be the set of haplotypes with fraction > 0.0 in the recursion so far.
                        // // For haplotype h, let similar_L(h) provide the set of haplotypes h' in L with dist(h,h') < x.
                        // // For haplotype h, let similar_R(h) provide the set of haplotypes h' in R with dist(h,h') < x.
                        // // A haplotype h may be explored with fraction > 0.0 if |similar_L(h)| >= |similar_R(h)|.
                        // // However, similar_*(h) should never return h itself.

                        // //achtung: the distance matrix should be gotten rid of the distance 0 entries, otherwise this will not work.
                        // //todo: double check that this works as expected.

                        // //find the current haplotype
                        // let current_haplotype = &self.haplotypes[haplotype_index];
                        // dbg!(&current_haplotype);

                        // //loop over haplotypes in the distance matrix and find haplotypes that are at distance x to the current haplotype (similar_L(h))
                        // // if let Some(distance_matrix) = &self.distance_matrix {
                        // // dbg!(&distance_matrix);
                        // let similar_l = self
                        //     .similarl_map
                        //     .as_ref()
                        //     .unwrap()
                        //     .get(&current_haplotype)
                        //     .unwrap();
                        // dbg!(&similar_l);
                        // dbg!(&fractions);
                        // //loop over haplotypes in collected fractions and find haplotypes that are at distance x to the current haplotype (similar_R(current_haplotype))
                        // let mut similar_r = 0;
                        // let similar_haplotypes_set = self.similar_haplotypes.as_ref().unwrap().get(&current_haplotype).unwrap();
                        // dbg!(&similar_haplotypes_set);
                        // for (h, f) in self.haplotypes[0..haplotype_index]
                        //     .to_vec()
                        //     .iter()
                        //     .zip(fractions[0..haplotype_index].to_vec().iter())
                        // {
                        //     if similar_haplotypes_set.contains(&h) && f > &NotNan::new(0.0).unwrap() {
                        //         similar_r += 1;
                        //     }
                        //     //for each collected fraction, loop over distance matrix and check if there is a haplotype with fraction > 0.0
                        //     // for ((h1, h2), distance) in distance_matrix.iter() {
                        //     //     if ((h1 == h && h2 == current_haplotype)
                        //     //         || (h2 == h && h1 == current_haplotype))
                        //     //         && (*distance < self.distance_threshold.unwrap())
                        //     //         && (f > &NotNan::new(0.0).unwrap())
                        //     //     {
                        //     //         // dbg!(&h, &f, &h1, &h2);
                        //     //         similar_r += 1;
                        //     //     }
                        //     // }
                        // }
                        // dbg!(&similar_r);
                        // //block the path in case similar_l is similar_R.
                        // //This way the recursion will continue only if similar_l >= similar_r.
                        // if *similar_l < similar_r {
                        //     dbg!(&"path is blocked");
                        //     return LogProb::ln_zero();
                        // }
                        // // }
                    // }
                }
                // dbg!(&fractions);
                let mut fractions = fractions.to_vec();
                fractions.push(fraction);
                self.calc_marginal(data, haplotype_index + 1, &mut fractions, joint_prob)
            };

            if haplotype_index == self.n_haplotypes - 1 {
                // let second_if = density(fraction_upper_bound);
                // second_if //clippy gives a warning
                density(fraction_upper_bound)
            } else {
                if fraction_upper_bound == NotNan::new(0.0).unwrap() {
                    // let last_if = density(NotNan::new(0.0).unwrap());
                    // last_if
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
                    } else if genotypes[i] == VariantStatus::Unknown || covered[i as u64] {
                        ()
                    }
                    // else if covered[i as u64] {
                    //     ()
                    // }
                    else if genotypes[i] == VariantStatus::NotPresent && !covered[i as u64] {
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
                if *fraction == NotNan::new(0.0).unwrap()
                    || *fraction == NotNan::new(0.5).unwrap()
                    || *fraction == NotNan::new(1.0).unwrap()
                {
                    prior_prob += LogProb::from(Prob(1.0 / 3.0))
                }
                // else if *fraction == NotNan::new(0.5).unwrap() {
                //     prior_prob += LogProb::from(Prob(1.0 / 3.0))
                // } else if *fraction == NotNan::new(1.0).unwrap() {
                //     prior_prob += LogProb::from(Prob(1.0 / 3.0))
                // }
                else {
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
