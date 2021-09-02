use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use failure::{format_err, Error};
use flate2::read::MultiGzDecoder;

fn open_gzip(filename: impl AsRef<Path>) -> Result<MultiGzDecoder<BufReader<File>>, Error> {
    let f = File::open(filename)?;
    Ok(MultiGzDecoder::new(BufReader::new(f)))
}

pub(crate) fn open_maybe_gz(filename: impl AsRef<Path>) -> Result<Box<dyn Read>, Error> {
    if filename.as_ref().exists() {
        open_by_extension(filename)
    } else {
        let mut filename_gz: PathBuf = filename.as_ref().into();
        let mut extension: OsString = filename_gz
            .extension()
            .map(|x| x.to_os_string())
            .unwrap_or("".into());
        extension.push(".gz");
        filename_gz.set_extension(extension);
        println!("opening: {:?}", filename_gz);

        if filename_gz.exists() {
            open_by_extension(filename_gz)
        } else {
            Err(format_err!(
                "File not found. '{:?}' and '{:?}' don't exist",
                filename.as_ref(),
                filename_gz
            ))
        }
    }
}

fn open_by_extension(filename: impl AsRef<Path>) -> Result<Box<dyn Read>, Error> {
    let res = match filename
        .as_ref()
        .extension()
        .and_then(std::ffi::OsStr::to_str)
    {
        Some("gz") => Box::new(open_gzip(filename)?) as Box<dyn Read>,
        _ => Box::new(File::open(filename)?) as Box<dyn Read>,
    };

    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn no_gzip() {
        let _ = open_maybe_gz("test/test_no_gzip.txt").expect("couldn't open file");
    }

    #[test]
    fn with_gzip() {
        let _ = open_maybe_gz("test/test_gzip.txt").expect("couldn't open file");
    }
}
