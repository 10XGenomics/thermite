use std::cmp::Ordering;
use std::str;
use std::str::FromStr;

use nom::{
    self,
    branch::alt,
    bytes::complete::{is_not, tag, take_until, take_while, take_while1},
    character::complete::{char, digit1},
    combinator::{all_consuming, map_res, opt},
    error::{ErrorKind, ParseError},
    sequence::{delimited, pair, separated_pair, terminated, tuple},
    IResult,
};

use failure::{format_err, Error};
use smallvec::SmallVec;

/// A parsed GTF record. The string
/// fields are borrowed from a buffer
/// containing the GTF line.
#[derive(Debug, Clone, PartialEq)]
pub struct Record<'a> {
    pub seqname: &'a [u8],
    pub source: &'a [u8],
    pub feature_type: &'a [u8],
    pub start: u64,
    pub end: u64,
    pub score: Option<f64>,
    pub strand: &'a [u8],
    pub frame: &'a [u8],
    pub attributes: AttrVec<'a>,
}

impl<'a> Record<'a> {
    pub fn all_attributes(&self) -> Vec<(String, String)> {
        let mut vec: Vec<(String, String)> = Vec::new();

        for (k, v) in &self.attributes {
            let k = std::str::from_utf8(k).unwrap().to_string();
            let v = std::str::from_utf8(v).unwrap().to_string();
            vec.push((k, v));
        }

        vec
    }

    pub fn get_attr(&self, attribute: &str) -> Result<String, Error> {
        for (k, v) in &self.attributes {
            if k == &attribute.as_bytes() {
                return Ok(std::str::from_utf8(v)?.to_string());
            }
        }

        Err(format_err!("attribute not found: {}", attribute))
    }
}

/// Parse one line of a GTF file into a `Record<'a>`. The
/// record will borrow slices from the input line.
pub fn parse_gtf_line<'a>(line: &'a [u8]) -> IResult<&[u8], Record<'a>> {
    let gtf_line = tuple((
        is_not("\t\r\n"), //chrom
        char('\t'),
        is_not("\t\r\n"), // source
        char('\t'),
        is_not("\t\r\n"), //chrom
        char('\t'),
        parse_u64, //start
        char('\t'),
        parse_u64, // end
        char('\t'),
        score, // score
        char('\t'),
        is_not("\t\r\n"), //strand -- FIXME!!
        char('\t'),
        is_not("\t\r\n"), //frame -- FIXME 0,1,2,or '.'
        char('\t'),
        gtf_attributes, //attributes
    ));

    let v = map_res(gtf_line, convert_to_record);
    all_consuming(v)(line)
}

impl<'a> PartialOrd for Record<'a> {
    fn partial_cmp(&self, other: &Record) -> Option<Ordering> {
        let r = self.seqname.cmp(&other.seqname);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.start.cmp(&other.start);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.end.cmp(&other.end);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.feature_type.cmp(&other.feature_type);
        if r != Ordering::Equal {
            return Some(r);
        }

        return Some(Ordering::Equal);
    }
}

/// convert ascii byte slice contain a decimal integer to u64
fn u64_from_str(input: &[u8]) -> Result<u64, std::num::ParseIntError> {
    u64::from_str(str::from_utf8(input).unwrap())
}

/// parse an integer from the input
fn parse_u64(input: &[u8]) -> IResult<&[u8], u64> {
    map_res(digit1, u64_from_str)(input)
}

/// return None unconditionally
fn empty_f64_option(_: &[u8]) -> Result<Option<f64>, std::num::ParseFloatError> {
    Ok(None)
}

/// convert ascii byte from
fn f64_option(input: &[u8]) -> Result<Option<f64>, std::num::ParseFloatError> {
    f64::from_str(str::from_utf8(input).unwrap()).map(|x| Some(x))
}

/// parse a GTF score field to Option<64>. An empty field '.',
/// will return None, otherwise Some(f64).
fn score(input: &[u8]) -> IResult<&[u8], Option<f64>> {
    alt((
        map_res(tag("."), empty_f64_option),
        map_res(nom::number::complete::recognize_float, f64_option),
    ))(input)
}

/// Is character a valid token for a GTF
/// attribute
#[inline]
fn is_token(c: u8) -> bool {
    match c {
        128..=255 => false,
        0..=31 => false,
        b'(' => false,
        b')' => false,
        b'<' => false,
        b'>' => false,
        b'@' => false,
        b',' => false,
        b';' => false,
        b':' => false,
        b'\\' => false,
        b'"' => false,
        b'/' => false,
        b'[' => false,
        b']' => false,
        b'?' => false,
        b'=' => false,
        b'{' => false,
        b'}' => false,
        b' ' => false,
        _ => true,
    }
}

type AttrVec<'a> = SmallVec<[(&'a [u8], &'a [u8]); 16]>;

fn gtf_attributes<'a>(input: &'a [u8]) -> IResult<&'a [u8], AttrVec<'a>> {
    terminated(
        separated_list_smallvec(
            pair(tag(";"), take_while1(|c| c == b' ')),
            separated_pair(
                take_while1(is_token),
                take_while1(|c| c == b' '),
                alt((
                    delimited(char('"'), take_until("\""), char('"')),
                    take_while1(is_token),
                )),
            ),
        ),
        opt(pair(opt(tag(";")), take_while(|c| c == b' '))),
    )(input)
}

/// raw fields of a GTF line, separated by tab characters
/// this is used transiently and will be converted to a
/// `Record<'a>`.
type RecInnerSep<'a> = (
    &'a [u8],
    char,
    &'a [u8],
    char,
    &'a [u8],
    char,
    u64,
    char,
    u64,
    char,
    Option<f64>,
    char,
    &'a [u8],
    char,
    &'a [u8],
    char,
    AttrVec<'a>,
);

fn convert_to_record<'a>(inp: RecInnerSep<'a>) -> Result<Record<'a>, std::num::ParseFloatError> {
    Ok(Record {
        seqname: inp.0,
        source: inp.2,
        feature_type: inp.4,
        start: inp.6,
        end: inp.8,
        score: inp.10,
        strand: inp.12,
        frame: inp.14,
        attributes: inp.16,
    })
}

/// Replacement for `separated_list` in nom, that returns items in a `SamllVec`
/// to avoid allocations in the tight inner loop of GTF attribute parsing.
fn separated_list_smallvec<I, O, O2, E, F, G>(
    sep: G,
    f: F,
) -> impl Fn(I) -> IResult<I, SmallVec<[O; 16]>, E>
where
    I: Clone + PartialEq,
    F: Fn(I) -> IResult<I, O, E>,
    G: Fn(I) -> IResult<I, O2, E>,
    E: ParseError<I>,
{
    use nom::Err;

    move |mut i: I| {
        let mut res = SmallVec::new();

        match f(i.clone()) {
            Err(Err::Error(_)) => return Ok((i, res)),
            Err(e) => return Err(e),
            Ok((i1, o)) => {
                if i1 == i {
                    return Err(Err::Error(E::from_error_kind(i1, ErrorKind::SeparatedList)));
                }

                res.push(o);
                i = i1;
            }
        }

        loop {
            match sep(i.clone()) {
                Err(Err::Error(_)) => return Ok((i, res)),
                Err(e) => return Err(e),
                Ok((i1, _)) => {
                    if i1 == i {
                        return Err(Err::Error(E::from_error_kind(i1, ErrorKind::SeparatedList)));
                    }

                    match f(i1.clone()) {
                        Err(Err::Error(_)) => return Ok((i, res)),
                        Err(e) => return Err(e),
                        Ok((i2, o)) => {
                            if i2 == i {
                                return Err(Err::Error(E::from_error_kind(
                                    i2,
                                    ErrorKind::SeparatedList,
                                )));
                            }

                            res.push(o);
                            i = i2;
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const ORIG: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const NO_QUOTES: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const BAD1: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene    _id ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const TWO_SPACE: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id  ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const GRCH_120: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2";"#;
    const TRAILING_SEMI_SPACE: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2"; "#;
    const TRAILING_SPACE: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2" "#;
    const EMPTY_ATTR_VAL: &[u8] = br#"chr21	HAVANA	gene	34073578	34106260	.	+	.	gene_id ""; gene_type "protein_coding"; gene_name ""; hgnc_id "HGNC:11038"; tag "overlapping_locus"; tag "retrogene"; havana_gene "OTTHUMG00000065821.4";"#;

    #[test]
    fn test_gtf_spaces() {
        let original = parse_gtf_line(ORIG);
        let no_quotes = parse_gtf_line(NO_QUOTES);
        assert_eq!(original, no_quotes);

        let two_space = parse_gtf_line(TWO_SPACE);
        assert_eq!(original, two_space);
    }

    #[test]
    fn bad_in_attr_name() {
        let b = parse_gtf_line(BAD1);
        println!("{:?}", b);
        assert!(b.is_err());
    }

    #[test]
    fn trailing_semi() {
        let b = parse_gtf_line(GRCH_120);
        println!("{:?}", b);
        assert!(b.is_ok());
    }

    #[test]
    fn trailing_space() {
        let b = parse_gtf_line(TRAILING_SPACE);
        println!("{:?}", b);
        assert!(b.is_ok());
    }

    #[test]
    fn trailing_semi_space() {
        let b = parse_gtf_line(TRAILING_SEMI_SPACE);
        println!("{:?}", b);
        assert!(b.is_ok());
    }

    #[test]
    fn empty_attr_value() {
        let b = parse_gtf_line(EMPTY_ATTR_VAL);
        println!("{:?}", b);
        assert!(b.is_ok());
    }

    /*
    // use for testing random GTF files

    use std::path::Path;
    use failu
    re::{Error, format_err};
    #[test]
    fn chr21_gtf() {
        parse_gtf_file("/Users/patrick/chr_21_genes.gtf").unwrap();
    }

    fn parse_gtf_file(gtf_path: impl AsRef<Path>) -> Result<(), Error> {
        use std::io::BufRead;

        let in_gtf = std::io::BufReader::new(std::fs::File::open(gtf_path)?);
        for (line_num, line) in in_gtf.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }

            let _rec = match crate::reference::parse_gtf::parse_gtf_line(&line.as_bytes()) {
                Ok((_, rec)) => rec,
                Err(e) => {
                    let msg = format_err!("Error parsing GTF on line: {}, {:?}", line_num, e);
                    return Err(msg);
                }
            };
        }

        Ok(())
    }
    */
}
