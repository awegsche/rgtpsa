use num_complex::ComplexFloat;
use num_traits::{Float, Num, One, Zero};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub};

use crate::{tpsa::TPSA, tpsa_utils::get_n_coeffs};

// -------------------------------------------------------------------------------------------------
// ---- Operators with Self ------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// ---- Add ----------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
impl<const NV: usize, const MO: usize, N: ComplexFloat> Add<&TPSA<NV, MO, N>> for &TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn add(self, rhs: &TPSA<NV, MO, N>) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> Add<TPSA<NV, MO, N>> for TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn add(self, rhs: TPSA<NV, MO, N>) -> Self::Output {
        let mut new_tpsa = self;
        new_tpsa += rhs;
        new_tpsa
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> AddAssign<TPSA<NV, MO, N>>
    for TPSA<NV, MO, N>
{
    fn add_assign(&mut self, rhs: TPSA<NV, MO, N>) {
        *self += &rhs;
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> AddAssign<&TPSA<NV, MO, N>>
    for TPSA<NV, MO, N>
{
    fn add_assign(&mut self, rhs: &TPSA<NV, MO, N>) {
        unsafe {
            for i in 0..self.n_coeffs() {
                *self.coeffs.get_unchecked_mut(i) =
                    *self.coeffs.get_unchecked_mut(i) + *rhs.coeffs.get_unchecked(i);
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
// ---- Sub ----------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
impl<const NV: usize, const MO: usize, N: ComplexFloat> Sub<TPSA<NV, MO, N>> for TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn sub(self, rhs: TPSA<NV, MO, N>) -> Self::Output {
        let mut coeffs = Vec::with_capacity(self.n_coeffs());

        unsafe {
            for i in 0..self.n_coeffs() {
                coeffs.push(*self.coeffs.get_unchecked(i) - *rhs.coeffs.get_unchecked(i));
            }
        }

        Self { coeffs }
    }
}

// -------------------------------------------------------------------------------------------------
// ---- Mul ----------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn mul_inplace(a: &Self, b: &Self, result: &mut Self) {
        unsafe {
            for &(i, j, k) in Self::get_mul_map() {
                *result.coeffs.get_unchecked_mut(k as usize) = *result.coeffs.get_unchecked(k as usize)
                    + *a.coeffs.get_unchecked(i as usize) * *b.coeffs.get_unchecked(j as usize);
            }
        }
    }
}
impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<&TPSA<NV, MO, N>> for &TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn mul(self, rhs: &TPSA<NV, MO, N>) -> Self::Output {
        let n_coeffs = get_n_coeffs(NV, MO);
        let mut result = Self::Output {
            coeffs: vec![N::zero(); n_coeffs],
        };

        TPSA::mul_inplace(self, rhs, &mut result);

        result
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<TPSA<NV, MO, N>> for &TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn mul(self, rhs: TPSA<NV, MO, N>) -> Self::Output {
        self * &rhs
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<&TPSA<NV, MO, N>> for TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn mul(self, rhs: &TPSA<NV, MO, N>) -> Self::Output {
        &self * rhs
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<TPSA<NV, MO, N>> for TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn mul(self, rhs: TPSA<NV, MO, N>) -> Self::Output {
        self * &rhs
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> MulAssign<&TPSA<NV, MO, N>>
    for TPSA<NV, MO, N>
{
    fn mul_assign(&mut self, rhs: &TPSA<NV, MO, N>) {
        *self = &*self * rhs;
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn powi(&self, n: i64) -> Self {
        // TODO: can be sped up by the squaring algorithm
        let mut result = Self::new(&[N::one()]);

        for _ in 0..n {
            result *= self;
        }

        result
    }
}

// -------------------------------------------------------------------------------------------------
// ---- Operators with N ---------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
impl<const NV: usize, const MO: usize, N: ComplexFloat> DivAssign<N> for TPSA<NV, MO, N> {
    fn div_assign(&mut self, rhs: N) {
        for a in self.coeffs.iter_mut() {
            *a = *a / rhs;
        }
    }
}
impl<const NV: usize, const MO: usize, N: ComplexFloat> Div<N> for TPSA<NV, MO, N> {
    type Output = Self;

    fn div(self, rhs: N) -> Self::Output {
        let mut new_tpsa = self;
        new_tpsa /= rhs;
        new_tpsa
    }
}

impl<const NV: usize, const MO: usize, N: ComplexFloat> MulAssign<N> for TPSA<NV, MO, N> {
    fn mul_assign(&mut self, rhs: N) {
        for a in self.coeffs.iter_mut() {
            *a = *a * rhs;
        }
    }
}
impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<N> for TPSA<NV, MO, N> {
    type Output = Self;

    fn mul(self, rhs: N) -> Self::Output {
        let mut new_tpsa = self;
        new_tpsa *= rhs;
        new_tpsa
    }
}
impl<const NV: usize, const MO: usize, N: ComplexFloat> Mul<N> for &TPSA<NV, MO, N> {
    type Output = TPSA<NV, MO, N>;

    fn mul(self, rhs: N) -> Self::Output {
        let new_tpsa = self.clone();
        new_tpsa * rhs
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn tpsa_1plus1() {
        type T = TPSA<3, 2, f64>;
        let t1 = T::new(&[1.0]);
        let t2 = T::new(&[1.0]);

        assert_eq!(T::new(&[2.0]), t1 + t2);
    }

    #[test]
    fn tpsa_add() {
        type T = TPSA<3, 2, f64>;
        let t1 = T::new(&[1.0, 0.0, 0.0, 1.0]);
        let t2 = T::new(&[1.0, 0.0, 1.0, 0.0]);
        let ts = T::new(&[2.0, 0.0, 1.0, 1.0]);

        assert_eq!(ts, t1 + t2);
    }

    #[test]
    fn tpsa_multiply() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy
        let t1 = T::new(&[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let t2 = T::new(&[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);

        assert_eq!(ts, t1 * t2);
    }

    #[test]
    fn tpsa_multiply_2x_3x() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy
        let t1 = T::new(&[0.0, 2.0, 0.0, 0.0, 0.0, 0.0]);
        let t2 = T::new(&[0.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 0.0, 0.0, 0.0, 6.0, 0.0]);

        assert_eq!(ts, t1 * t2);
    }

    #[test]
    fn tpsa_x_square() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy
        let t1 = T::new(&[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);

        assert_eq!(ts, t1.powi(2));
    }

    #[test]
    fn tpsa_binomial() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy  y^2
        let t1 = T::new(&[0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0]);

        assert_eq!(ts, t1.powi(2));
    }

    #[test]
    fn tpsa_half() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy  y^2
        let t1 = T::new(&[0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0]);

        assert_eq!(ts, t1 / 2.0);
    }

    #[test]
    fn tpsa_double() {
        type T = TPSA<3, 2, f64>;
        //                  1    x    y    z  x^2   xy  y^2
        let t1 = T::new(&[0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let ts = T::new(&[0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

        assert_eq!(ts, t1 * 2.0);
    }
    use rand::prelude::*;
}
