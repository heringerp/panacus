pub struct WalkSplitter<'a> {
    slice: &'a str,
}

impl<'a> WalkSplitter<'a> {
    #[inline]
    pub fn new(slice: &'a str) -> Self {
        Self { slice }
    }
}

impl<'a> Iterator for WalkSplitter<'a> {
    type Item = &'a str;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.slice.is_empty() {
            return None;
        }

        let next_marker_idx = self.slice[1..]
            .chars()
            .position(|b| b == '>' || b == '<')
            .map(|i| i + 1)
            .unwrap_or(self.slice.len());

        let (token, rest) = self.slice.split_at(next_marker_idx);
        self.slice = rest;
        Some(token)
    }
}

pub struct WalkByteSplitter<'a> {
    slice: &'a [u8],
}

impl<'a> WalkByteSplitter<'a> {
    #[inline]
    pub fn new(slice: &'a [u8]) -> Self {
        Self { slice }
    }
}

impl<'a> Iterator for WalkByteSplitter<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.slice.is_empty() {
            return None;
        }

        let next_marker_idx = self.slice[1..]
            .iter()
            .position(|&b| b == b'>' || b == b'<')
            .map(|i| i + 1)
            .unwrap_or(self.slice.len());

        let (token, rest) = self.slice.split_at(next_marker_idx);
        self.slice = rest;
        Some(token)
    }
}
