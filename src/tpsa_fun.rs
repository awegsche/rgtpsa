use crate::tpsa::TPSA;
use num_complex::ComplexFloat;
use num_traits::One;
use num_traits::Zero;

impl<const NV: usize, const MO: usize, N: ComplexFloat> TPSA<NV, MO, N> {
    pub fn sin(&self) -> TPSA<NV, MO, N> {
        let mut result = self.clone();
        let mut factor = self.clone();
        let self_squared = self * self;

        for k in 1..=((MO as i32) + 2) {
            let factorial = -(2.0 * k as f64) * (2.0 * k as f64 + 1.0);
            factor *= &self_squared;
            //factor *= N::from(-1.0).unwrap();
            factor /= N::from(factorial).unwrap();
            result += &factor;
        }
        result
    }

    pub fn cos(&self) -> TPSA<NV, MO, N> {
        let mut result = Self::one();
        let mut factor = Self::one();

        for k in 1..=((MO as u32) * 2) {
            let factorial = N::from((2 * k - 1) * (2 * k)).unwrap();
            factor *= self;
            factor *= self;
            factor *= N::from(-1.0).unwrap();
            factor /= factorial;
            result += &factor;
        }
        result
    }

    pub fn exp(&self) -> TPSA<NV, MO, N> {
        let mut result = Self::one();
        let mut factor = Self::one();

        for k in 1..=((MO as u32) * 2) {
            let factorial = N::from(k).unwrap();
            factor *= self;
            factor /= factorial;
            result += &factor;
        }
        result
    }

    /// Derivation w.r.t. the specified variable
    ///
    /// _Note_: the variable index is in human-readable format (i.e. starting at 1).
    ///
    /// # Example
    ///
    /// ```
    /// # use gtpsa::tpsa::TPSA;
    /// # type T = TPSA<3,2,f64>;
    /// // t = 2x
    /// let t = T::from_x(2.0);
    ///
    /// // a = t^2 = 4x^2
    /// let a = t.powi(2);
    ///
    /// // da/dx = 8x
    /// let da_dx = a.der(1);
    /// assert_eq!(da_dx, T::from_x(8.0));
    /// // d^2a/dx^2 = 8
    /// let d2a_dx2 = da_dx.der(1);
    /// assert_eq!(d2a_dx2, T::from_1(8.0));
    /// ```
    pub fn der(&self, ivar: usize) -> TPSA<NV, MO, N> {
        let mut result = Self::zero();

        let map = Self::get_der_map();
        let slice = (NV - ivar) * self.n_coeffs();

        for (i, x) in self.coeffs.iter().enumerate() {
            let (idx, factor) = map[slice + i];
            if idx == usize::MAX || factor.is_zero() {
                continue;
            }
            result.coeffs[idx] = *x * N::from(factor).unwrap();
        }

        result
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn derive_x_squared() {
        type T = TPSA<3, 2, f64>;
        let t = T::from_x(3.0);
        let t = t.powi(2);

        let dt_dx = t.der(1);
        assert_eq!(dt_dx, T::from_x(18.0));

        let d2t_dx2 = dt_dx.der(1);
        assert_eq!(d2t_dx2, T::from_1(18.0));
    }
}
