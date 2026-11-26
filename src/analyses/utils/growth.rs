use rayon::prelude::*;

use crate::{file_formats::gfa_parser::ThresholdContainer, hist::Hist, util::Threshold};

/// Predicts continued growth using Smoothed Good-Toulmin (Orlitsky et al. 2016)
/// Target is a float multiplying the original sample size to the new sample size (including all original samples). So target = 2.0 means we predict the growth up to a sample size of 2.0 * n
/// The returned Vec goes from n+1 ... n * target + 1
#[expect(dead_code, reason = "Not yet implemented growth curve extrapolation")]
pub fn predict_ul(hist: &Hist, target: f64, max_rel_se: Option<f64>) -> Vec<f64> {
    let max_rel_se = max_rel_se.unwrap_or(0.01);
    let n = hist.get_max_frequency() as f64;
    let m = 1..((target - 1.0) * n + 1.0) as usize;
    let t = m.map(|x| x as f64 / n);
    let growths = t.map(|x| predict_sgt(hist, x, max_rel_se));
    growths.collect()
}

/// Calculates just the last predicted growth value instead of the full array. Use this everytime
/// you don't need the full array (e.g. regional-growth) as lots of calculations can be skipped
/// then.
pub fn predict_ul_at_t(hist: &Hist, target: f64, max_rel_se: Option<f64>) -> f64 {
    let max_rel_se = max_rel_se.unwrap_or(0.01);
    let t = target - 1.0;
    predict_sgt(hist, t, max_rel_se).floor()
}

fn sgt_at_int_l(hist: &Hist, t: f64, l: usize, c_buf: &mut Vec<f64>) -> (f64, f64) {
    if t <= 0.0 || l == 0 || hist.get_max_frequency() < 1 {
        return (0.0, 0.0);
    }

    if t <= 1.0 {
        let (u_hat, var_hat, _) = hist.get_hist_values().iter().skip(1).enumerate().fold(
            (0.0_f64, 0.0_f64, 1.0_f64),
            |(u, var, w_prev), (idx, &phi)| {
                // w = t^i (since this is carried over from each iteration to the next)
                let w = w_prev * t;
                let phi = phi as f64;
                let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
                (u + sign * w * phi, var + w * w * phi, w)
            },
        );
        if u_hat < 0.0 {
            return (0.0, var_hat);
        }
        return (u_hat, var_hat);
    }

    // Forward pass: c[j - 1] = nCr(l, j) * 2^j - (1 - q)^l
    let one_minus_q = t / (t + 2.0);
    let mut c_curr = one_minus_q.powi(l as i32);
    c_buf.clear();
    // Reserves additional space (but since it was cleared before it only increases by
    // l - current_capacity)
    c_buf.reserve(l);
    for j in 1..=l {
        c_curr *= 2.0 * (l - j + 1) as f64 / j as f64;
        c_buf.push(c_curr);
    }

    // Backward pass: w_i = c_i + w_{i+1} / t
    let mut w = 0.0_f64;
    let mut u_hat = 0.0_f64;
    let mut var_hat = 0.0_f64;
    for idx in (0..l).rev() {
        w = c_buf[idx] + w / t;
        if idx < hist.get_max_frequency() {
            let phi = hist[idx + 1] as f64;
            let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
            u_hat += sign * w * phi;
            var_hat += w * w * phi;
        }
    }
    if u_hat < 0.0 {
        u_hat = 0.0;
    }
    (u_hat, var_hat)
}

fn predict_sgt(hist: &Hist, t: f64, max_rel_se: f64) -> f64 {
    let s_obs: f64 = hist
        .get_hist_values()
        .iter()
        .skip(1)
        .map(|&x| x as f64)
        .sum();
    let mut c_buf = Vec::new();

    if t <= 1.0 {
        let (u_hat, _) = sgt_at_int_l(hist, t, 1, &mut c_buf);
        return s_obs + u_hat;
    }

    let u_inc: f64 = hist
        .get_hist_values()
        .iter()
        .enumerate()
        .map(|(idx, &x)| idx as f64 * x as f64)
        .sum();
    let l_float = 1.0_f64.max(0.5 * (u_inc * t * t / (t - 1.0)).log(3.0));
    let l_low = l_float.floor() as usize;
    let alpha = l_float - (l_low as f64);

    let (mut best_u, mut best_var) = sgt_at_int_l(hist, t, 1, &mut c_buf);
    for l in 2..=l_low {
        let (u_curr, var_curr) = sgt_at_int_l(hist, t, l, &mut c_buf);
        let rel_se = var_curr.sqrt() / 1.0_f64.max(s_obs + u_curr);
        if !rel_se.is_finite() || rel_se > max_rel_se {
            return s_obs + best_u;
        }
        (best_u, best_var) = (u_curr, var_curr);
    }

    let (u_high, var_high) = sgt_at_int_l(hist, t, l_low + 1, &mut c_buf);
    let u_blend = (1.0 - alpha) * best_u + alpha * u_high;
    let var_blend = (1.0 - alpha).powi(2) * best_var + alpha.powi(2) * var_high;

    let blend_rel_se = var_blend.sqrt() / 1.0_f64.max(s_obs + u_blend);
    if blend_rel_se.is_finite() && blend_rel_se <= max_rel_se {
        return s_obs + u_blend;
    }
    s_obs + best_u
}

pub fn calc_growth(hist: &Hist, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage();

    if n > 0 {
        let quorum = usize::max(1, t_quorum.to_absolute(n));
        if quorum == 1 {
            calc_growth_union(hist, t_coverage)
        } else if quorum >= n {
            calc_growth_core(hist, t_coverage)
        } else {
            calc_growth_quorum(hist, t_coverage, t_quorum)
        }
    } else {
        Vec::new()
    }
}

pub fn calc_all_growths(
    hist: &Hist,
    hist_aux: &ThresholdContainer,
    insert_zero: bool,
) -> Vec<Vec<f64>> {
    let mut growths: Vec<Vec<f64>> = hist_aux
        .coverage
        .par_iter()
        .zip(&hist_aux.quorum)
        .map(|(c, q)| {
            log::info!(
                "calculating growth for coverage >= {} and quorum >= {}",
                &c,
                &q
            );
            calc_growth(hist, c, q)
        })
        .collect();
    // insert empty row for 0 element
    if insert_zero {
        for g in &mut growths {
            g.insert(0, f64::NAN);
        }
    }
    growths
}

fn calc_growth_union(hist: &Hist, t_coverage: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from 0..n (both included)
    let c = usize::max(1, t_coverage.to_absolute(n));

    let mut pangrowth: Vec<f64> = vec![0.0; n];
    let mut n_fall_m: f64 = 0.0;
    let tot = hist.get_hist_values()[c..].iter().sum::<usize>() as f64;

    // perc_mult[i] contains the percentage of combinations that
    // have an item of multiplicity i
    let mut perc_mult: Vec<f64> = Vec::with_capacity(n + 1);
    perc_mult.resize(n + 1, 0.0);

    for m in 1..n + 1 {
        let mut y: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for (i, perc_mult_value) in perc_mult.iter_mut().enumerate().take(n - m + 1).skip(c) {
            *perc_mult_value += (n as f64 - m as f64 - i as f64 + 1.0).log2();
            y += ((hist.get_hist_values()[i] as f64).log2() + *perc_mult_value - n_fall_m).exp2();
        }

        pangrowth[m - 1] = tot - y;
    }

    pangrowth
}

fn calc_growth_core(hist: &Hist, t_coverage: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from 0..n (both included)
    let c = usize::max(1, t_coverage.to_absolute(n + 1));
    let mut n_fall_m: f64 = 0.0;
    let mut pangrowth: Vec<f64> = vec![0.0; n];

    // In perc_mult[i] is contained the percentage of combinations
    // that have an item of multiplicity i
    let mut perc_mult: Vec<f64> = Vec::with_capacity(n + 1);
    perc_mult.resize(n + 1, 0.0);

    for m in 1..n + 1 {
        let mut y: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for (i, perc_mult_value) in perc_mult
            .iter_mut()
            .enumerate()
            .take(n + 1)
            .skip(usize::max(m, c))
        {
            *perc_mult_value += (i as f64 - m as f64 + 1.0).log2();
            y += ((hist.get_hist_values()[i] as f64).log2() + *perc_mult_value - n_fall_m).exp2();
        }
        pangrowth[m - 1] = y;
    }

    pangrowth
}

fn calc_growth_quorum(hist: &Hist, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from [0..n]
    let c = usize::max(1, t_coverage.to_absolute(n));
    let quorum = t_quorum.to_relative(n);
    let mut pangrowth: Vec<f64> = vec![0.0; n];

    let mut n_fall_m: f64 = 0.0;
    let mut m_fact: f64 = 0.0;

    let mut perc_mult: Vec<f64> = vec![0.0; n + 1];
    let mut q: Vec<Vec<f64>> = vec![vec![0.0; n + 1]; n + 1];

    for m in 1..n + 1 {
        m_fact += (m as f64).log2();
        let m_quorum = (m as f64 * quorum).ceil() as usize;

        //100% quorum
        let mut yl: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for (i, perc_mult_value) in perc_mult
            .iter_mut()
            .enumerate()
            .take(n + 1)
            .skip(usize::max(m, c))
        {
            *perc_mult_value += (i as f64 - m as f64 + 1.0).log2();
            yl += ((hist.get_hist_values()[i] as f64).log2() + *perc_mult_value - n_fall_m).exp2();
        }

        //[m_quorum, 100) quorum
        let mut yr: f64 = 0.0;
        for (i, row) in q.iter_mut().enumerate().take(n).skip(m_quorum) {
            let mut sum_q = 0.0;
            let mut add = false;
            for (j, value) in row
                .iter_mut()
                .enumerate()
                .take(m)
                .skip(usize::max(m_quorum, c))
            {
                if n + j + 1 > i + m && j <= i {
                    if *value == 0.0 {
                        *value = choose(i, j);
                    }
                    *value += (n as f64 - i as f64 - m as f64 + 1.0 + j as f64).log2();
                    *value -= (m as f64 - j as f64).log2();
                    sum_q += (*value + m_fact - n_fall_m).exp2();
                    add = true;
                }
            }
            if add {
                yr += ((hist.get_hist_values()[i] as f64).log2() + sum_q.log2()).exp2();
            }
        }
        pangrowth[m - 1] = yl + yr;
    }
    pangrowth
}

fn choose(n: usize, k: usize) -> f64 {
    let mut res: f64 = 0.0;
    if k > n {
        return 0.0;
    }

    let k = if k > n - k { n - k } else { k };

    let n = n as f64;

    for i in 0..k {
        res += (n - i as f64).log2();
        res -= (i as f64 + 1.0).log2();
    }
    res
}

#[cfg(test)]
mod tests {

    use crate::util::Threshold;

    use super::*;

    #[test]
    fn test_sgt() {
        let hist = Hist::from_values(vec![0, 33, 16, 14, 11, 8, 6]);
        let growths: Vec<_> = predict_ul(&hist, 3.0, None)
            .iter()
            .map(|x| x.floor() as i32)
            .collect();
        assert_eq!(
            growths,
            vec![93, 97, 101, 105, 108, 110, 112, 114, 116, 118, 119, 121]
        );
    }

    #[test]
    fn test_hist_calc_growth_union() {
        let hist = Hist::from_values(vec![0, 5, 3, 2]);

        let t_coverage = Threshold::Absolute(0);
        let test_growth: Vec<f64> = vec![5.666666666666667, 8.333333333333334, 10.0];
        let growth = calc_growth_union(&hist, &t_coverage);
        assert_eq!(growth, test_growth, "Wrong growth union");
    }

    #[test]
    fn test_hist_calc_growth_core() {
        let hist = Hist::from_values(vec![0, 5, 3, 2]);

        let t_coverage = Threshold::Absolute(0);
        let test_core: Vec<f64> = vec![5.666666666666666, 3.0, 2.0];
        let core = calc_growth_core(&hist, &t_coverage);
        assert_eq!(core, test_core, "Wrong growth core");
    }

    #[test]
    fn test_hist_calc_growth_quorum() {
        let hist = Hist::from_values(vec![0, 5, 3, 2, 3, 5, 0, 4, 2, 1]);

        let t_coverage = Threshold::Absolute(0);
        let t_quorum = Threshold::Relative(0.9);
        let test_growth: Vec<f64> = vec![
            11.88888888888889,
            7.027777777777777,
            4.761904761904761,
            3.4444444444444438,
            2.5952380952380953,
            2.0,
            1.5555555555555545,
            1.2222222222222217,
            1.0,
        ];
        let growth = calc_growth_quorum(&hist, &t_coverage, &t_quorum);
        assert_eq!(growth, test_growth, "Wrong growth quorum");
    }
}
