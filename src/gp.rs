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

impl FromStr for Direction {
    type Err = ParseGraphPosError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => { Ok(Direction::Forward) },
            "-" => { Ok(Direction::RevComp) },
            _   => { Err(ParseGraphPosError) }
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
        let (input, gp) = alt((Self::parse_graphpos, Self::parse_graphpos_incomplete))(s)
            .map_err(|_| ParseGraphPosError)?;
        if !input.is_empty() { return Err(ParseGraphPosError); }
        Ok(gp)
    }

}

use nom::branch::alt;
use nom::character::complete::digit1;
use nom::character::complete::one_of;
use nom::bytes::complete::tag;
use nom::IResult;
impl GraphPos {
    fn parse_graphpos(input: &str) -> IResult<&str, Self> {
        let (input, id) = digit1(input)?;
        let (input, sign) = one_of("+-")(input)?;
        let (input, _) = tag(":")(input)?;
        let (input, pos) = digit1(input)?;
        let gp = GraphPos{
            // these unwraps are safe because of nom
            id   : id.parse().unwrap(),
            sign : sign.to_string().parse().unwrap(),
            pos  : pos.parse().unwrap()
        };
        return Ok((input, gp));
    }

    fn parse_graphpos_incomplete(input: &str) -> IResult<&str, Self> {
        let (input, id) = digit1(input)?;
        let (input, sign) = one_of("+-")(input)?;
        let gp = GraphPos{
            id : id.parse().unwrap(),
            sign: sign.to_string().parse().unwrap(),
            pos: 0,
        };
        return Ok((input, gp));
    }
}

#[test]
fn can_parse_graph_pos() {
    let gp: GraphPos = "10+".parse().unwrap();
    assert_eq!(gp, GraphPos{id: 10, sign: Direction::Forward, pos: 0});
    let gp: GraphPos = "10+:5".parse().unwrap();
    assert_eq!(gp, GraphPos{id: 10, sign: Direction::Forward, pos: 5});
    let gp = "10+:".parse::<GraphPos>();
    assert!(gp.is_err());
}

impl GraphPos {
    pub fn to_path(self) -> String {
        let mut repr = String::new();
        match self.sign {
            Direction::Forward => { repr.push('>'); },
            Direction::RevComp => { repr.push('<'); }
        }
        write!(&mut repr, "{}", self.id).unwrap();
        return repr;
    }
}
