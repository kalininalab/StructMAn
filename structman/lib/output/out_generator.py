class OutputGenerator:
    null_symbol = '-'

    def __init__(self):
        self.columns = []
        self.current_line = []
        self.header_map = {}

    def add_headers(self, headers):
        n = len(self.columns)
        for pos, header in enumerate(headers):
            self.header_map[header] = pos + n
            self.columns.append(header)

    def get_header(self):
        header = '\t'.join(self.columns) + '\n'
        return header

    def add_value(self, header, value):
        if value is None:
            value = OutputGenerator.null_symbol
        if not isinstance(value, str):
            value = str(value)

        pos = self.header_map[header]
        if len(self.current_line) <= pos:
            self.current_line += [OutputGenerator.null_symbol] * (1 + (pos - len(self.current_line)))
        self.current_line[pos] = value

    def pop_line(self):
        if len(self.current_line) < len(self.columns):
            self.current_line += [OutputGenerator.null_symbol] * (len(self.columns) - len(self.current_line))
        line = '\t'.join(self.current_line) + '\n'
        self.current_line = []
        return line
