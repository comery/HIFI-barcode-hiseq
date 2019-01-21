class Error(Exception):
    pass

class Line(str):
    """A line of text with associated filename and line number."""
    def error(self, message):
        """Return an error relating to this line."""
        return Error("{0}({1}): {2}\n{3}"
                     .format(self.filename, self.lineno, message, self))

class Lines(object):
    """Lines(filename, iterator) wraps 'iterator' so that it yields Line
    objects, with line numbers starting from 1. 'filename' is used in
    error messages.

    """
    def __init__(self, filename, iterator):
        self.filename = filename
        self.lines = enumerate(iterator, start=1)

    def __iter__(self):
        return self

    def __next__(self):
        lineno, s = next(self.lines)
        line = Line(s)
        line.filename = self.filename
        line.lineno = lineno
        return line

    # For compatibility with Python 2.
    next = __next__

def read_fastq(filename, iterator):
    """Read FASTQ data from 'iterator' (which may be a file object or any
    other iterator that yields strings) and generate tuples (sequence
    name, sequence data, quality data). 'filename' is used in error
    messages.

    """
    # This implementation follows the FASTQ specification given here:
    # <http://nar.oxfordjournals.org/content/38/6/1767.full>
    import re
    at_seqname_re = re.compile(r'@(.+)$')
    sequence_re = re.compile(r'[!-*,-~]*$')
    plus_seqname_re = re.compile(r'\+(.*)$')
    quality_re = re.compile(r'[!-~]*$')

    lines = Lines(filename, iterator)
    for line in lines:
        # First line of block is @<seqname>.
        m = at_seqname_re.match(line)
        if not m:
            raise line.error("Expected @<seqname> but found:")
        seqname = m.group(1)
        try:
            # One or more lines of sequence data.
            sequence = []
            for line in lines:
                m = sequence_re.match(line)
                if not m:
                    break
                sequence.append(m.group(0))
            if not sequence:
                raise line.error("Expected <sequence> but found:")

            # The line following the sequence data consists of a plus
            # sign and an optional sequence name (if supplied, it must
            # match the sequence name from the start of the block).
            m = plus_seqname_re.match(line)
            if not m:
                raise line.error("Expected +[<seqname>] but found:")
            if m.group(1) not in ['', seqname]:
                raise line.error("Expected +{} but found:".format(seqname))

            # One or more lines of quality data, containing the same
            # number of characters as the sequence data.
            quality = []
            n = sum(map(len, sequence))
            while n > 0:
                line = next(lines)
                m = quality_re.match(line)
                if not m:
                    raise line.error("Expected <quality> but found:")
                n -= len(m.group(0))
                if n < 0:
                    raise line.error("<quality> is longer than <sequence>:")
                quality.append(m.group(0))

            yield seqname, ''.join(sequence), ''.join(quality)

        except StopIteration:
            raise line.error("End of input before sequence was complete:")