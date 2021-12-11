from copy import copy, deepcopy
from itertools import chain


class Matrix(object):
    """

    Representation of matrices.

    """

    def __init__(self, entries):
        """

        Parameters:
        -----------
        entries: list
            Two dimensional list of entries.

        Notes:
        ------
        Entries have to be integers.

        """

        self.entries = entries
        self.num_rows = len(entries)
        self.num_cols = len(entries[0])

    def get_row(self, row_index):
        """

        Get a row from matrix.

        Parameters:
        -----------
        row_index: int
            Index of a row.

        Output:
        -------
        row: list
            A row.

        """

        return self.entries[row_index]

    def get_col(self, col_index):
        """

        Get a column from matrix.

        Parameters:
        -----------
        col_index: int
            Index of a column.

        Output:
        -------
        col: list
            A column.

        """

        col = []
        for row in self.entries:
            col.append(row[col_index])
        return col

    def change_row(self, row_index, values):
        """

        Return a matrix with changed values of row.

        Parameters:
        -----------
        matrix: list
            A matrix.
        row_index: int
            Index of a row.
        values: list
            A list with desired values of row.

        Output:
        -------
        matrix: list
            A matrix with changed rows.

        """

        self.entries[row_index] = values
        return self.entries

    def change_col(self, col_index, values):
        """

        Return a matrix with changed values of column.

        Parameters:
        -----------
        col_index: int
            Index of a column.
        values: list
            A list with desired values of column.

        Output:
        -------
        entries: list
            A entries of a matrix with changed values.

        """

        for row_index in range(self.num_rows):
            self.entries[row_index][col_index] = values[col_index]
        return self.entries

    def add_rows(self, index_row_to_change, index_row_to_add):
        """

        Add one row to another.

        Parameters:
        -----------
        index_row_to_change: int
            Index of changed row.
        index_row_to_add: int
            Index of added row to another.

        Notes:
        ------
        The method uses modulo 2 addition of entries.

        """

        row_to_change = self.get_row(index_row_to_change)
        row_to_add = self.get_row(index_row_to_add)
        summed_row = []
        for i in range(self.num_cols):
            summed_entries = (row_to_change[i] + row_to_add[i]) % 2
            summed_row.append(summed_entries)
        self.change_row(index_row_to_change, summed_row)

    def add_cols(self, index_col_to_change, index_col_to_add):
        """

        Add one col to another.

        Parameters:
        -----------
        index_col_to_change: int
            Index of changed col.
        index_col_to_add: int
            Index of added col to another.

        Notes:
        ------
        The method uses modulo 2 addition of entries.

        """

        col_to_change = self.get_col(index_col_to_change)
        col_to_add = self.get_col(index_col_to_add)
        summed_col = []
        for i in range(self.num_rows):
            summed_entries = (col_to_change[i] + col_to_add[i]) % 2
            summed_col.append(summed_entries)
        self.change_row(index_col_to_change, summed_col)

    def non_zero(self, entry_index):
        """

        Arranges rows and cols so that specified entry on diagonal is non zero.

        Parameters:
        -----------
        entry_index: int
            Index of a row (and column) that contains the entry.

        Output:
        -------
        value_check: bool
            Boolean that tells whether the entry is not zero.

        Notes:
        ------
        Intended to use only in smith_normal_form method due to algorithm
        optimizations.

        It will move a row, if non zero entry is on the same column as the
        desired place.

        It will move a columns, if non zero entry is right of
        the desired place.

        """

        zero_row = True
        for i in range(self.num_cols):
            if self.entries[entry_index][i] != 0:
                zero_row = False
        if zero_row is True:
            return False

        if self.entries[entry_index][entry_index] != 0:
            return True

        for row_index in chain(range(entry_index),
                               range(entry_index+1, self.num_rows)):

            if self.entries[row_index][entry_index] != 0:
                copied_row = copy(self.get_row(row_index))
                copied_entry_row = copy(self.get_row(entry_index))

                self.change_row(entry_index, copied_row)
                self.change_row(row_index, copied_entry_row)

                return True

        for col_index in chain(range(entry_index),
                               range(entry_index+1, self.num_cols)):
            if self.entries[entry_index][col_index] != 0:
                copied_col = copy(self.get_col(col_index))
                copied_entry_col = copy(self.get_col(entry_index))

                self.change_col(entry_index, copied_col)
                self.change_col(col_index, copied_entry_col)

                return True

        return False

    def snf(self):
        """

        Calculate Smith normal form of a matrix.

        Output:
        -------
        snf: matrix
            Smith normal form of a matrix.

        Notes:
        ------
        The method uses modulo 2 addition.

        The method is not consistent with theory, but it is good enough
        for getting rank of a matrix.

        """

        snf = deepcopy(self)
        max_dim = min(snf.num_rows, snf.num_cols)

        for dim in range(max_dim):
            value_check = snf.non_zero(dim)

            if value_check is True:
                for row_index in chain(range(dim), range(dim+1, snf.num_rows)):
                    if snf.entries[row_index][dim] == 1:
                        snf.add_rows(row_index, dim)

        return snf

    def snf_rank(self):
        """

        Calculate rank of Smith normal form of a matrix.

        Output:
        -------
        snf_rank: int
            Rank of Smith normal form of a matrix.

        Notes:
        ------
        The method uses modulo 2 addition, with exception of calculating rank
        itself.

        """

        snf = self.snf()
        max_dim = min(snf.num_rows, snf.num_cols)

        snf_rank = 0
        for i in range(max_dim):
            snf_rank += snf.entries[i][i]

        return snf_rank
