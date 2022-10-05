use num_complex::ComplexFloat;
use std::{
    collections::HashMap,
    fmt::Debug,
    fmt::Display,
    ops::{Index, IndexMut},
};

use num_traits::{One, Zero};
use once_cell::sync::OnceCell;

use crate::tpsa_utils::{get_exponents_for_order_flat, get_n_coeffs};

// -------------------------------------------------------------------------------------------------
// ---- TPSA struct definition ---------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

/// TPSA struct
#[derive(PartialEq, Eq, Clone)]
pub struct TPSA<const NV: usize, const MO: usize, N: ComplexFloat> {
    pub(crate) coeffs: Vec<N>,
}

// -------------------------------------------------------------------------------------------------
// ---- Initialisation of struct and const statics -------------------------------------------------
// -------------------------------------------------------------------------------------------------

impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn new(coeffs: &[N]) -> Self {
        let n_coeffs = get_n_coeffs(NV, MO);
        assert!(n_coeffs >= coeffs.len());

        let mut a = vec![N::zero(); n_coeffs];

        a[..coeffs.len()].clone_from_slice(coeffs);

        Self { coeffs: a }
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn n_coeffs(&self) -> usize {
        self.coeffs.len()
    }

    pub(crate) fn get_exp() -> &'static Vec<u8> {
        static INSTANCE: OnceCell<Vec<u8>> = OnceCell::new();

        INSTANCE.get_or_init(|| {
            let mut stack = Vec::with_capacity(get_n_coeffs(NV, MO));
            for order in 0..=MO {
                get_exponents_for_order_flat([0; NV], 0, order, &mut stack);
            }
            stack
        })
    }

    pub(crate) fn get_index_from_exponents_map() -> &'static HashMap<&'static [u8], usize> {
        static INSTANCE: OnceCell<HashMap<&'static [u8], usize>> = OnceCell::new();

        INSTANCE.get_or_init(|| {
            let mut temp_map = HashMap::new();

            let exps = Self::get_exp();
            let n_coeffs = get_n_coeffs(NV, MO);

            for idx in 0..n_coeffs {
                temp_map.insert(&exps[idx * NV..(idx + 1) * NV], idx);
            }
            temp_map
        })
    }

    /// returns the index range into the map from [`get_exp`](get_exp)
    /// # Example
    /// ```
    /// #
    /// # struct TPSA {};
    /// # impl TPSA {
    /// #   pub(crate) fn get_exp_idx_range(idx: usize) -> std::ops::Range<usize> {
    /// #       (idx * NV)..(idx+1)*NV
    /// #   }
    /// #   pub(crate) fn get_exp() -> Vec<u8> {
    /// #       vec![0;10]
    /// #   }
    /// # }
    /// #
    /// # const NV: usize = 3;
    /// #
    /// impl TPSA {
    ///     fn priv_function() {
    ///         let exp = Self::get_exp();
    ///         let exponent: &[u8] = &exp[Self::get_exp_idx_range(0)];
    ///     }
    /// }
    /// ```
    #[inline]
    pub(crate) fn get_exp_idx_range(idx: usize) -> std::ops::Range<usize> {
        (idx * NV)..(idx + 1) * NV
    }
    #[inline]
    pub fn get_mul_map() -> &'static Vec<(u16, u16, u16)> {
        static INSTANCE: OnceCell<Vec<(u16, u16, u16)>> = OnceCell::new();

        INSTANCE.get_or_init(|| {
            let temp_map = Self::get_index_from_exponents_map();
            let n_coeffs = get_n_coeffs(NV, MO);
            let exps = Self::get_exp();

            let mut map = Vec::new();

            assert!(n_coeffs < u16::MAX.into());

            for idx1 in 0..n_coeffs {
                let exp1 = &exps[Self::get_exp_idx_range(idx1)];
                for idx2 in 0..n_coeffs {
                    let exp2 = &exps[Self::get_exp_idx_range(idx2)];
                    let mut sum = [0; NV];
                    let mut ord_sum = 0;
                    for i in 0..NV {
                        sum[i] = exp1[i] + exp2[i];
                        ord_sum += sum[i] as usize;
                    }
                    if ord_sum <= MO {
                        map.push((idx1 as u16, idx2 as u16, temp_map[&sum[..]] as u16));
                    }
                }
            }

            map
        })
    }

    #[inline]
    pub(crate) fn get_der_map() -> &'static Vec<(usize, f64)> {
        static INSTANCE: OnceCell<Vec<(usize, f64)>> = OnceCell::new();

        INSTANCE.get_or_init(|| {
            let exp_map = Self::get_index_from_exponents_map();
            let n_coeffs = get_n_coeffs(NV, MO);
            let exps = Self::get_exp();

            let mut map = Vec::new();

            for var in 0..NV {
                for i in 0..n_coeffs {
                    let exponents = &exps[Self::get_exp_idx_range(i)];
                    let mut der_exponents = [0; NV];
                    der_exponents.copy_from_slice(exponents);
                    if der_exponents[var] > 0 {
                        let factor = der_exponents[var] as f64;
                        der_exponents[var] -= 1;
                        let idx = exp_map[&der_exponents[..]];
                        map.push((idx, factor));
                    } else {
                        map.push((usize::MAX, 0.0));
                    }
                }
            }

            for i in 0..10 {
                dbg!(map[i]);
            }
            map
        })
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> One for TPSA<NV, MO, N> {
    fn one() -> Self {
        Self::new(&[N::one()])
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> Zero for TPSA<NV, MO, N> {
    fn zero() -> Self {
        Self::new(&[N::zero()])
    }

    fn is_zero(&self) -> bool {
        for a in self.coeffs.iter() {
            if !a.is_zero() {
                return false;
            }
        }
        true
    }
}

// -------------------------------------------------------------------------------------------------
// ---- Printing -----------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
impl<const NV: usize, const MO: usize, N: ComplexFloat + std::fmt::LowerExp> TPSA<NV, MO, N> {
    pub(crate) fn print(
        &self,
        print_zeros: bool,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        let exponents = Self::get_exp();
        write!(f, "{:>3} | {:25} | EXP", "I", "COEFFICIENT")?;
        for _ in 0..(NV - 1).max(0) {
            write!(f, "   ")?;
        }
        writeln!(f, "| ORDER")?;

        write!(f, "----+---------------------------+-")?;
        for _ in 0..NV {
            write!(f, "---")?;
        }
        writeln!(f, "+------")?;

        for idx in 0..self.n_coeffs() {
            let a = self.coeffs[idx];
            if !print_zeros && a.is_zero() {
                continue;
            }
            write!(f, "{:3} | {:25.18e} | ", idx, self.coeffs[idx])?;
            let mut order = 0;
            for e in 1..=NV {
                let exp = exponents[(idx + 1) * NV - e];
                order += exp;
                write!(f, "{:2} ", exp)?;
            }

            writeln!(f, "| {:2}", order)?;
        }
        Ok(())
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat + std::fmt::LowerExp> Display
    for TPSA<NV, MO, N>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.print(false, f)
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat + std::fmt::LowerExp> Debug
    for TPSA<NV, MO, N>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.print(true, f)
    }
}
// -------------------------------------------------------------------------------------------------
// ---- Accessors ----------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//

impl<const NV: usize, const MO: usize, N: ComplexFloat> Index<usize> for TPSA<NV, MO, N> {
    type Output = N;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coeffs[index]
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> IndexMut<usize> for TPSA<NV, MO, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coeffs[index]
    }
}

// -------------------------------------------------------------------------------------------------
// ---- Convenience --------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn from_1(a: N) -> Self {
        Self::new(&[a])
    }

    /// Creates a TPSA with the first variable filled with `a`
    ///
    /// Example:
    /// ```
    /// # use rgtpsa::tpsa::TPSA;
    /// # type T = TPSA<3,2,f64>;
    /// let x = T::from_x(1.0);
    ///
    /// //  this should be equal to 0 + 1*x:
    /// assert_eq!(x, T::new(&[0.0, 1.0]));
    /// ```
    pub fn from_x(a: N) -> Self {
        Self::new(&[N::zero(), a])
    }

    pub fn from_y(a: N) -> Self {
        Self::new(&[N::zero(), N::zero(), a])
    }

    pub fn from_z(a: N) -> Self {
        Self::new(&[N::zero(), N::zero(), N::zero(), a])
    }

    pub fn from_var_n(n: usize, a: N) -> Self {
        let mut result = Self::zero();
        result[n] = a;
        result
    }
}

#[cfg(test)]
mod test {
    use super::*;
}
