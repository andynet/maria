use std::fmt::Display;
use std::num::ParseIntError;
use std::str::FromStr;

#[derive(Clone, Copy, Eq, PartialEq, Hash)]
pub enum Direction {
    Forward,
    RevComp
}

use std::fmt;
impl Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Direction::Forward => { write!(f, "+") },
            Direction::RevComp => { write!(f, "-") }
        }
    } 
}

#[derive(Clone, Eq, Hash, PartialEq)]
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
        let direction;
        match s.bytes().last().ok_or(ParseGraphPosError)? {
            b'+' => { direction = Direction::Forward; },
            b'-' => { direction = Direction::RevComp; },
            _    => { return Err(ParseGraphPosError); }
        }
        let node_pos = 0;
        Ok(GraphPos { id: node_id, sign: direction, pos: node_pos })
    }
}

