use std::fmt::Display;
use std::fmt::Write;
use std::fmt;
use std::num::ParseIntError;
use std::str::FromStr;

#[derive(Clone, Copy, Eq, PartialEq, Hash, Debug, Default)]
pub enum Direction {
    #[default]
    Forward,
    RevComp
}

impl Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Direction::Forward => { write!(f, "+") },
            Direction::RevComp => { write!(f, "-") }
        }
    } 
}

#[derive(Clone, Eq, Hash, PartialEq, Debug, Default, Copy)]
pub struct GraphPos {
    pub id: usize,
    pub sign: Direction,
    pub pos: usize,
}

#[derive(Debug)]
pub struct ParseGraphPosError;
impl From<ParseIntError> for ParseGraphPosError {
    fn from(_: ParseIntError) -> Self { ParseGraphPosError }
}

impl FromStr for GraphPos {
    type Err = ParseGraphPosError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let node_id: usize = s[..s.len()-1].parse()?;
        let direction = match s.bytes().last().ok_or(ParseGraphPosError)? {
            b'+' => { Direction::Forward },
            b'-' => { Direction::RevComp },
            _    => { return Err(ParseGraphPosError); }
        };
        let node_pos = 0;
        Ok(GraphPos { id: node_id, sign: direction, pos: node_pos })
    }
}

impl GraphPos {
    pub fn to_path(self) -> String {
        let mut repr = String::new();
        match self.sign {
            Direction::Forward => { repr.push('>'); },
            Direction::RevComp => { repr.push('<'); }
        }
        // repr.write_fmt(format_args!("{}", self.id));
        write!(&mut repr, "{}", self.id).unwrap();
        return repr;
    }
}
