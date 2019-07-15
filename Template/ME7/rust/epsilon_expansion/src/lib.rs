use bit_array::BitArray;
use generic_array::{ArrayLength, GenericArray};
use generic_array_legacy::ArrayLength as ArrayLengthLegacy;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use typenum::{Integer, NonZero, Unsigned, B1, U32};

#[macro_use]
extern crate derivative;

#[derive(Derivative)]
#[derivative(Debug, Clone(bound = ""))]
/// An ε-expansion with `U` elements, where `M` is the highest pole in ε.
/// For example:
/// ```
/// use epsilon_expansion::EpsilonExpansion;
/// EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[(-3, 1.0), (0, 2.0), (2, 3.0)]);
/// ```
/// creates an ε-expansion with range `1/ε^6` to `ε^6` and reads `1/ε^3+2+3ε^2`.
pub struct EpsilonExpansion<U: Unsigned + NonZero, M: Integer>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    coeffs: GenericArray<f64, U>,
    nonzero: BitArray<u32, U>,
    _phantom: std::marker::PhantomData<M>,
}

impl<U: Unsigned + NonZero, M: Integer> PartialEq for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

impl<U: Unsigned + NonZero, M: Integer> EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    #[inline]
    /// Create an empty ε-expansion.
    pub fn new() -> EpsilonExpansion<U, M> {
        EpsilonExpansion {
            coeffs: GenericArray::default(),
            nonzero: BitArray::new(),
            _phantom: std::marker::PhantomData,
        }
    }

    #[inline]
    /// Create an ε-expansion from a slice where each entry in `elem`
    /// is a tuple of a power in ε and a coefficient.
    pub fn from_slice(elem: &[(isize, f64)]) -> EpsilonExpansion<U, M> {
        let mut e = EpsilonExpansion::new();
        for (pow, coeff) in elem {
            e.set(*pow, *coeff);
        }
        e
    }

    #[inline]
    /// Set the `ε^pow`th coefficient to `coeff`. This function
    /// checks if `pow` is in range.
    pub fn set(&mut self, pow: isize, coeff: f64) {
        if pow < M::to_isize() {
            eprintln!("truncating {}/ε^{}", coeff, -pow);
            return;
        }
        let i = (pow - M::to_isize()) as usize;
        if i >= self.coeffs.len() {
            eprintln!("truncating {}*ε^{}", coeff, pow);
            return;
        }

        self.nonzero.set(i, coeff != 0.);
        self.coeffs[i] = coeff;
    }

    #[inline]
    /// Add `coeff` to the `ε^pow`th coefficient This function
    /// checks if `pow` is in range.
    pub fn add_coeff(&mut self, pow: isize, coeff: f64) {
        if pow < M::to_isize() {
            eprintln!("truncating {}/ε^{}", coeff, -pow);
            return;
        }
        let i = (pow - M::to_isize()) as usize;
        if i >= self.coeffs.len() {
            eprintln!("truncating {}*ε^{}", coeff, pow);
            return;
        }

        self.coeffs[i] += coeff;
        self.nonzero.set(i, self.coeffs[i] != 0.);
    }
}

impl<U: Unsigned + NonZero, M: Integer> fmt::Display for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, (oc, nz)) in self.coeffs.iter().zip(&self.nonzero).enumerate() {
            if nz {
                write!(f, "+")?;
                oc.fmt(f)?;
                if (i as isize) < -M::to_isize() {
                    write!(f, "/ε^{}", -M::to_isize() - i as isize)?;
                } else if i as isize > -M::to_isize() {
                    write!(f, "ε^{}", i as isize + M::to_isize())?;
                }
            }
        }
        write!(f, "")
    }
}

impl<U: Unsigned + NonZero, M: Integer> fmt::LowerExp for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, (oc, nz)) in self.coeffs.iter().zip(&self.nonzero).enumerate() {
            if nz {
                write!(f, "+")?;
                fmt::LowerExp::fmt(oc, f)?;
                if (i as isize) < -M::to_isize() {
                    write!(f, "/ε^{}", -M::to_isize() - i as isize)?;
                } else if i as isize > -M::to_isize() {
                    write!(f, "ε^{}", i as isize + M::to_isize())?;
                }
            }
        }
        write!(f, "")
    }
}

impl<U: Unsigned + NonZero, M: Integer> Neg for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn neg(self) -> EpsilonExpansion<U, M> {
        let mut res = (*self).clone();

        for (oc, nz) in res.coeffs.iter_mut().zip(&self.nonzero) {
            if nz {
                *oc = (*oc).neg();
            }
        }
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Neg for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn neg(mut self) -> EpsilonExpansion<U, M> {
        for (oc, nz) in self.coeffs.iter_mut().zip(&self.nonzero) {
            if nz {
                *oc = (*oc).neg();
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<&EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        let mut res = (*self).clone();
        res.nonzero.union(&other.nonzero); // TODO: update after 0-check?

        for (rc, (oc, nz)) in res
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc += oc;
            }
        }
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(self, mut other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        other.nonzero.union(&self.nonzero);

        for (rc, (oc, nz)) in other
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(&self.nonzero))
        {
            if nz {
                *rc += oc;
            }
        }
        other
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<&EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(mut self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        self.nonzero.union(&other.nonzero);

        for (rc, (oc, nz)) in self
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc += oc;
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(mut self, other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        self.nonzero.union(&other.nonzero);

        for (rc, (oc, nz)) in self
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc += oc;
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<f64> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(mut self, other: f64) -> EpsilonExpansion<U, M> {
        self.add_coeff(0, other);
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Add<f64> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn add(self, other: f64) -> EpsilonExpansion<U, M> {
        let mut res = self.clone();
        res.add_coeff(0, other);
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<&EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        let mut res = (*self).clone();
        res.nonzero.union(&other.nonzero);

        for (rc, (oc, nz)) in res
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc -= oc;
            }
        }
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(self, mut other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        other.nonzero.union(&self.nonzero);

        for (rc, (oc, nz)) in other
            .coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(&self.nonzero))
        {
            if nz {
                *rc = oc - *rc;
            } else {
                *rc = -*rc;
            }
        }
        other
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<&EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(mut self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        self.nonzero.union(&other.nonzero);

        for (rc, (oc, nz)) in self
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc -= oc;
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(mut self, other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        self.nonzero.union(&other.nonzero);

        for (rc, (oc, nz)) in self
            .coeffs
            .iter_mut()
            .zip(other.coeffs.iter().zip(&other.nonzero))
        {
            if nz {
                *rc -= oc;
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<f64> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(mut self, other: f64) -> EpsilonExpansion<U, M> {
        self.add_coeff(0, -other);
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Sub<f64> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn sub(self, other: f64) -> EpsilonExpansion<U, M> {
        let mut res = self.clone();
        res.add_coeff(0, -other);
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<&EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        let mut res = EpsilonExpansion::new();

        for (i, (oc, nz)) in self.coeffs.iter().zip(&self.nonzero).enumerate() {
            if nz {
                for (j, (oc1, nz1)) in other.coeffs.iter().zip(&other.nonzero).enumerate() {
                    if nz1 {
                        res.add_coeff(
                            i as isize + M::to_isize() + j as isize + M::to_isize(),
                            oc * oc1,
                        );
                    }
                }
            }
        }
        res
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<EpsilonExpansion<U, M>> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(self, other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        self.mul(&other)
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<&EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(self, other: &EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        (&self).mul(other)
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<EpsilonExpansion<U, M>> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(self, other: EpsilonExpansion<U, M>) -> EpsilonExpansion<U, M> {
        (&self).mul(&other)
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<f64> for EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(mut self, other: f64) -> EpsilonExpansion<U, M> {
        for (oc, nz) in self.coeffs.iter_mut().zip(&self.nonzero) {
            if nz {
                *oc *= other;
            }
        }
        self
    }
}

impl<U: Unsigned + NonZero, M: Integer> Mul<f64> for &EpsilonExpansion<U, M>
where
    U: ArrayLength<f64>,
    U: Add<U32>,
    <U as Add<U32>>::Output: Sub<B1>,
    <<U as Add<U32>>::Output as Sub<B1>>::Output: Div<U32>,
    <<<U as Add<U32>>::Output as Sub<B1>>::Output as Div<U32>>::Output: ArrayLengthLegacy<u32>,
{
    type Output = EpsilonExpansion<U, M>;
    fn mul(self, other: f64) -> EpsilonExpansion<U, M> {
        let mut res = self.clone();
        for (oc, nz) in res.coeffs.iter_mut().zip(&res.nonzero) {
            if nz {
                *oc *= other;
            }
        }
        res
    }
}

#[test]
fn test() {
    use EpsilonExpansion;

    let a =
        EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[(-3, 1.0), (0, 2.0), (2, 3.0)]);
    let b =
        EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[(-1, 2.0), (1, 5.0), (2, 6.0)]);
    let c = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-3, 1.0),
        (-1, 2.0),
        (0, 2.0),
        (1, 5.0),
        (2, 9.0),
    ]);
    let d = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-4, 2.0),
        (-2, 5.0),
        (-1, 10.0),
        (1, 16.0),
        (2, 12.0),
        (3, 15.),
        (4, 18.),
    ]);
    let e =
        EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[(-3, 1.0), (0, 8.0), (2, 3.0)]);
    let f = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-3, 6.0),
        (0, 12.0),
        (2, 18.0),
    ]);
    let h = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-3, 1.0),
        (-1, -2.0),
        (0, 2.0),
        (1, -5.0),
        (2, -3.0),
    ]);
    let i = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-3, 1.0),
        (0, -4.0),
        (2, 3.0),
    ]);
    let j = EpsilonExpansion::<typenum::U13, typenum::N6>::from_slice(&[
        (-3, -1.0),
        (0, -2.0),
        (2, -3.0),
    ]);

    assert_eq!(&a + &b, c);
    assert_eq!(&a + b.clone(), c);
    assert_eq!(a.clone() + &b, c);
    assert_eq!(a.clone() + b.clone(), c);

    assert_eq!(&a - &b, h);
    assert_eq!(&a - b.clone(), h);
    assert_eq!(a.clone() - &b, h);
    assert_eq!(a.clone() - b.clone(), h);

    assert_eq!(&a * &b, d);
    assert_eq!(&a * b.clone(), d);
    assert_eq!(a.clone() * &b, d);
    assert_eq!(a.clone() * b.clone(), d);

    assert_eq!(&a + 6.0, e);
    assert_eq!(a.clone() + 6.0, e);

    assert_eq!(&a - 6.0, i);
    assert_eq!(a.clone() - 6.0, i);

    assert_eq!(&a * 6.0, f);
    assert_eq!(a.clone() * 6.0, f);

    assert_eq!(-&a, j);
    assert_eq!(-a.clone(), j);
}
