import re
from itertools import islice

class LazyLineObject:
    """
    Instead of storing ``lines`` as an array, this object can be used.
    It reduces the memory usage drastically! It looks up lines only when needed.
    """
    def __init__(self, file, start, end):
       self.file = file
       self.start = start
       self.end = end

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f"LazyLineObject for file {self.file}, lines {self.start}-{self.end}"

    def __repr__(self):
        return f"LazyLineObject for file {self.file}, lines {self.start}-{self.end}"

    def __iter__(self):
        with open(self.file, "r") as lines:
            for line in islice(lines, self.start, self.end + 1):
                yield line.rstrip("\n")

    def __getitem__(self, key):
        if key >= len(self):
            raise KeyError("key too big")
        with open(self.file, "r") as lines:
            for line in islice(lines, self.start + key, self.start + key + 1):
                return line.rstrip()

    def full_text(self):
        text = ""
        with open(self.file, "r") as lines:
            for line in islice(lines, self.start, self.end + 1):
                text += line.rstrip() + "\n"
        return text

    def search_for_block(self, start, end, count=1, join=" ", max_len=1000, format_line=None):
        """
        Search through a file (lines) and locate a block starting with "start" (inclusive) and ending with "end" (exclusive).

        Args:
            start (str): a pattern that matches the start of the block (can contain special characters)
            end (str): a pattern that matches the end of the block (can contain special characters) - ``None`` removes this (so a selection of ``max_lines`` is guaranteed)
            count (int): how many matches to search for
            join (str): spacer between lines
            max_len (int): maximum length of matches (to prevent overflow)
            format_line (function): function to perform to each line before adding to match (e.g. remove leading space)

        Returns:
            a single match (str) if count == 1 or a list of matches (str) if count > 1.
        """
        assert isinstance(count, int), "count needs to be an integer"
        assert isinstance(max_len, int), "count needs to be an integer"
        assert isinstance(join, str), "join needs to be a string"

        if count == 0:
            return None

        current_match = ""
        current_len = 0
        match = [None] * count

        #### we want a regex that will never match anything - and quickly - so trying to match something before the start of the line works
        if end is None:
            end = "a^"

        start_pattern = re.compile(start)
        end_pattern = re.compile(end)

        index = 0
        for line in self:
            if current_match:
                if end_pattern.search(line) or current_len >= max_len:
                    match[index] = current_match
                    current_match = None
                    index += 1
                    current_len = 0

                    if index == count:
                        break
                else:
                    if format_line is not None:
                        current_match = current_match + join + format_line(line.lstrip())
                    else:
                        current_match = current_match + join + line.lstrip()
                    current_len += 1
            else:
                if start_pattern.search(line):
                    if format_line is not None:
                        current_match = format_line(line.lstrip())
                    else:
                        current_match = line.lstrip()
                    current_len = 1

        if count == 1:
            return match[0]
        else:
            return match


    def find_parameter(self, parameter, expected_length, which_field, split_on=None, cast_to_float=True):
        """
        Args:
            parameter (string): test to search for
            expected_length (int): how many fields there should be
            which_field (int or list): which field(s) the parameter is (zero-indexed)
            split_on (str): additional non-space field on which to split
            cast_to_float (Bool): whether or not to cast extracted value to float
        Returns:
            a list of all the extracted values
        """
        if not isinstance(which_field, list):
            which_field = [which_field]

        if not isinstance(expected_length, int):
            raise TypeError("expected_length must be type int!")

        for n in which_field:
            if not isinstance(n, int):
                raise TypeError("which_field must be type int!")
            if n >= expected_length:
                raise ValueError("can't expect a field after the last field!")

        matches = []
        pattern = False

        try:
            pattern = re.compile(parameter)
        except Exception as e:
            raise ValueError("pattern {pattern} cannot be compiled as a regex; try again!")

        if pattern:
            for line in self:
                if pattern.search(line):
                    fields = re.split(" +", line)
                    if split_on:
                        fields2 = []
                        for field in fields:
                            fields2 = fields2 + field.split(split_on)
                        fields = fields2
                    fields = list(filter(None, fields))

                    if len(fields) == expected_length:
                        desired_fields = []
                        for n in which_field:
                            if cast_to_float:
                                try:
                                    desired_fields.append(float(fields[n]))
                                except:
                                    desired_fields.append(0)
                            else:
                                desired_fields.append(fields[n])
                        if len(desired_fields) == 1:
                            matches.append(desired_fields[0])
                        else:
                            matches.append(desired_fields)
            return matches

