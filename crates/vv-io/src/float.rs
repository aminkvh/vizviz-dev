//! Fast decimal parsing for the numbers structure files actually contain
//! (`-12.345`, `1.00`, `85`), with std as the fallback for anything else.

/// Parses an ASCII decimal without exponent. Falls back to `str::parse`
/// (exponents, `inf`) and returns `None` for non-numbers.
pub fn parse_f32(s: &[u8]) -> Option<f32> {
    let (negative, digits) = match s.first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    let mut int: u64 = 0;
    let mut i = 0;
    while i < digits.len() && digits[i].is_ascii_digit() {
        if i >= 15 {
            return fallback(s);
        }
        int = int * 10 + (digits[i] - b'0') as u64;
        i += 1;
    }
    let int_digits = i;
    let mut frac: u64 = 0;
    let mut scale: i32 = 0;
    if i < digits.len() && digits[i] == b'.' {
        i += 1;
        while i < digits.len() && digits[i].is_ascii_digit() {
            if scale < 9 {
                frac = frac * 10 + (digits[i] - b'0') as u64;
                scale += 1;
            }
            i += 1;
        }
    }
    if i != digits.len() || (int_digits == 0 && scale == 0) {
        return fallback(s);
    }
    let value = int as f64 + frac as f64 / 10f64.powi(scale);
    let value = value as f32;
    Some(if negative { -value } else { value })
}

pub fn parse_i32(s: &[u8]) -> Option<i32> {
    let (negative, digits) = match s.first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || digits.len() > 9 || !digits.iter().all(u8::is_ascii_digit) {
        return std::str::from_utf8(s).ok()?.trim().parse().ok();
    }
    let v = digits
        .iter()
        .fold(0i32, |acc, d| acc * 10 + (d - b'0') as i32);
    Some(if negative { -v } else { v })
}

fn fallback(s: &[u8]) -> Option<f32> {
    std::str::from_utf8(s).ok()?.trim().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn parses_typical_coordinates() {
        assert_eq!(parse_f32(b"17.047"), Some(17.047));
        assert_eq!(parse_f32(b"-0.500"), Some(-0.5));
        assert_eq!(parse_f32(b"1.00"), Some(1.0));
        assert_eq!(parse_f32(b"85"), Some(85.0));
        assert_eq!(parse_f32(b".5"), Some(0.5));
        assert_eq!(parse_f32(b"+3.25"), Some(3.25));
        assert_eq!(parse_f32(b"1e3"), Some(1000.0));
        assert_eq!(parse_f32(b"?"), None);
        assert_eq!(parse_f32(b""), None);
        assert_eq!(parse_f32(b"-"), None);
        assert_eq!(parse_i32(b"-12"), Some(-12));
        assert_eq!(parse_i32(b"1234567890"), Some(1234567890));
        assert_eq!(parse_i32(b"."), None);
    }

    proptest! {
        #[test]
        fn matches_std_within_one_ulp(int in -99999i32..99999, frac in 0u32..1000, scale in 0u32..4) {
            let text = if scale == 0 {
                format!("{int}")
            } else {
                format!("{int}.{:0width$}", frac % 10u32.pow(scale), width = scale as usize)
            };
            let expected: f32 = text.parse().unwrap();
            let got = parse_f32(text.as_bytes()).unwrap();
            let ulp = (expected.abs() * f32::EPSILON).max(f32::MIN_POSITIVE);
            prop_assert!((got - expected).abs() <= ulp, "{text}: {got} vs {expected}");
        }
    }
}
