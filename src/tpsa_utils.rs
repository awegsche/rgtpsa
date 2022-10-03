pub(crate) const fn factorial(n: usize) -> usize {
    let mut result = 1;
    let mut i = 2;
    while i <= n {
        result *= i;
        i+=1;
    }

    result
}
pub(crate) const fn get_n_coeffs(nv: usize, mo: usize) -> usize {
    let mut n = 1;

    let mut i = 1;
    while i <= (mo+nv) {
        n *= i;
        i+=1;
    }
    i = 1;
    while i <= mo {
        n /= i;
        i+=1;
    }
    i = 1;
    while i <= nv {
        n /= i;
        i+=1;
    }
    n
}

pub(crate) fn get_exponents_for_order_flat<const NV: usize>(
    last: [u8; NV],
    index: usize,
    sum: usize,
    stack: &mut Vec<u8>,
) {
    if sum == 0 {
        stack.extend_from_slice(&last);
        return;
    }
    if index == NV {
        //stack.extend_from_slice(&last);
        return;
    }
    for i in 0..=sum {
        let mut copylast = last;
        copylast[index] = i as u8;
        get_exponents_for_order_flat(copylast, index + 1, sum - i, stack);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(2), 2);
        assert_eq!(factorial(3), 6);
    }

    #[test]
    fn n_coeffs_3_0() {
        assert_eq!(get_n_coeffs(3, 0), 1);
    }

    #[test]
    fn n_coeffs_3_1() {
        assert_eq!(get_n_coeffs(3, 1), 4);
    }

    #[test]
    fn n_coeffs_3_2() {
        assert_eq!(get_n_coeffs(3, 2), 10);
    }

}
