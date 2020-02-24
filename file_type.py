import pandas as pd

class FileType:
    def __init__(self, file_name, file_extension):
        self.file_extension = file_extension
        self.file_name = file_name

    def read_file(self):
        obj = None
        full_file_name = self.file_name + self.file_extension

        if '.xlsx' == self.file_extension:
            # TODO replace sheet_name instead of Sheet1 name, and uncomment the line: sheet_name = input('What is the Sheet name in the excel?')
            # sheet_name = input('What is the Sheet name in the excel?')
            obj = pd.read_excel(full_file_name, 'Sheet1')
        elif '.csv' == self.file_extension:
            obj = pd.read_csv(self.file_name)
        elif '.txt' == self.file_extension:
            obj = pd.read_csv(self.file_name)

        return obj
