from itertools import combinations

import networkx as nx

from tda.matrix import Matrix
from tda.metrics import euclidean_metric
from tda.point import Point
from tda.utils import are_unique


class VietorisRipsComplex(object):
    """

    Representation of the Vietoris-Rips complex.

    Parameters:
    -----------
    points: list
        A list of `Point` objects.
    epsilon: float
        A positive real number.
    metric: callable
        A function that calculates distance between `Point` objects.

    References:
    -----------
    https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex

    """

    def __init__(self, points, epsilon, metric, validate_points=True):

        if epsilon <= 0:
            raise ValueError('Epsilon has to be greater than 0.')

        if not points:
            raise ValueError('List of points cannot be empty.')

        if validate_points:
            if not are_unique(points):
                raise ValueError('Points passed as input are not unique.')

        self.points = points
        self.epsilon = epsilon
        self.metric = metric
        self.n_points = len(points)

        self.graph = nx.Graph()
        self.n_edges = None
        self.faces = None
        self.simplices = None
        self.max_dim = None
        self.faces = None

    def create_graph(self):
        """

        Create a graph from points.

        """

        self._add_vertices()
        self._add_edges()

    def _add_vertices(self):
        """

        Add vertices to the graph.

        """

        self.graph.add_nodes_from(self.points)

    def _add_edges(self):
        """

        Add edges to the graph.

        Notes:
        ------
        The graph is constructed by placing
        an edge between each point such that
        the distance between them is less than
        epsilon.

        """

        self.n_edges = 0
        for i in range(len(self.points)):
            for j in range(i + 1, len(self.points)):
                p1, p2 = self.points[i], self.points[j]
                if self.metric(p1, p2) < self.epsilon:
                    self.graph.add_edge(p1, p2)
                    self.n_edges += 1

    def _remove_edges(self):
        """

        Remove edges from graph.

        """

        self.graph.clear_edges()

    def find_simplices(self):
        """

        Find simplices in graph.

        Output:
        -------
        simplices : tuple
            Tuple of maximal simplices in the complex.

        Notes:
        ------
        This function will not list simplices that are contained in other
        simplices.

        References:
        -----------
        https://en.wikipedia.org/wiki/Clique_(graph_theory)
        https://en.wikipedia.org/wiki/Simplex

        """

        self.simplices = tuple(map(tuple, nx.find_cliques(self.graph)))
        return self.simplices

    def find_faces(self):
        """

        Find faces.

        Output:
        -------
        faces : tuple
            Tuple of faces.

        References:
        -----------
        https://en.wikipedia.org/wiki/Clique_(graph_theory)
        https://en.wikipedia.org/wiki/Simplex

        """

        faces = set()
        for s in self.simplices:
            n_edges = len(s)
            for face_dim in range(n_edges, 0, -1):
                for face in combinations(s, face_dim):
                    curr_face = tuple(sorted(face, key=lambda x: x.name))
                    faces.add(curr_face)
        self.faces = tuple(faces)
        return self.faces

    def find_faces_with_dim(self, dim):
        """

        Find faces of a given dimension.

        Parameters:
        -----------
        dim : int
            A non-negative integer.

        Output:
        -------
        faces : tuple
            Tuple of faces of the given dimension.

        References:
        -----------
        https://en.wikipedia.org/wiki/Clique_(graph_theory)
        https://en.wikipedia.org/wiki/Simplex

        """

        if dim < 0:
            raise ValueError('A non-negative dimension was expected.')

        faces_with_dim = set()
        for s in self.simplices:
            if len(s) < dim:
                continue
            for face in combinations(s, dim + 1):
                curr_face = tuple(sorted(face, key=lambda x: x.name))
                faces_with_dim.add(curr_face)
        faces_with_dim = tuple(faces_with_dim)
        return faces_with_dim

    def change_epsilon(self, epsilon):
        """

        Change epsilon.

        Parameters:
        -----------
        epsilon: float
            A positive float.

        """

        if epsilon <= 0:
            raise ValueError('Epsilon has to be greater than 0.')

        self.epsilon = epsilon

        self.simplices = None
        self.n_edges = None
        self.faces = None
        self.simplices = None
        self.max_dim = None
        self.faces = None

        self._remove_edges()
        self._add_edges()

    @property
    def complex_dimension(self):
        """

        Calculate the dimension of the simplicial complex.

        Output:
        -------
        max_dim : int
            Dimension of the simplicial complex

        """

        self.max_dim = 0
        for s in self.simplices:
            self.max_dim = max(self.max_dim, len(s) - 1)
        return self.max_dim

    def check_nesting(self, higher_simplices, lower_simplices):
        """

        Find which simplices are nested in others.

        Parameters:
        -----------
        higher_simplices : tuple
            A tuple of simplices of higher dimension.
        lower_simplices : tuple
            A tuple of simplices of lower dimension.

        Output:
        -------
        nested_simplices : dict
            A dictionary of nested simplices.

        Notes:
        ------
        Simplices have to be of one higher dimension that the others.

        References:
        -----------
        https://en.wikipedia.org/wiki/Betti_number
        develop

        """

        nested = dict()
        for lower_simplex in lower_simplices:
            for higher_simplex in higher_simplices:
                is_nested = True
                for p in lower_simplex:
                    if p not in higher_simplex:
                        is_nested = False
                        break
                if is_nested:
                    if higher_simplex not in nested:
                        nested[higher_simplex] = list()
                    nested[higher_simplex].append(lower_simplex)
        return nested

    def boundary_operator_matrix(self, n):
        """

        Calculate matrix of a n-th boundary operator.

        Parameters:
        -----------
        n : int
            A non-negative integer.

        Output:
        -------
        boundary_matrix : matrix
            Matrix of a n-th boundary operator.
        matrix_rows: dict
            Dictionary of rows of the matrix.
        matrix_cols: dict
            Dictionary of columns of the matrix.

        Notes:
        ------
        The operation of chain group is addition with Z_2 coefficients.

        """

        if n < 0:
            raise ValueError('"n" has to be a non-negative integer.')

        higher_simplices = self.find_faces_with_dim(n+1)
        lower_simplices = self.find_faces_with_dim(n)

        boundary_dict = self.check_nesting(
            higher_simplices, lower_simplices)

        boundary_entries = []
        matrix_rows = {}
        reversed_matrix_rows = {}
        row_index = 0
        for lower_simplex in lower_simplices:
            boundary_entries.append([])
            matrix_rows[row_index] = lower_simplex
            reversed_matrix_rows[lower_simplex] = row_index
            row_index += 1

        matrix_cols = {}
        col_index = 0
        for higher_simplex in higher_simplices:
            for row in boundary_entries:
                row.append(0)

            list_of_faces = boundary_dict[higher_simplex]
            for face in list_of_faces:
                row_index = reversed_matrix_rows[face]
                boundary_entries[row_index][col_index] = 1

            matrix_cols[col_index] = higher_simplex
            col_index += 1

        boundary_matrix = Matrix(boundary_entries)

        return boundary_matrix, matrix_rows, matrix_cols

    def betti_numbers(self):
        """

        Calculate Betti numbers of the Vietoris-Rips complex.

        Outputs:
        --------
        betti_numbers : list
            List of Betti numbers of a complex.

        Notes:
        ------
        The operation of chain group is addition with Z_2 coefficients.

        References:
        -----------
        https://en.wikipedia.org/wiki/Betti_number
        https://youtu.be/gVq_xXnwV-4

        """

        raise NotImplementedError()

    def nth_betti_number(self, n):
        """

        Calculate n-th Betti number of the Vietoris-Rips complex.

        Outputs:
        --------
        nth-betti_number : int
            N-th Betti number of a complex.

        Notes:
        ------
        The operation of chain group is addition with Z_2 coefficients.

        References:
        -----------
        https://en.wikipedia.org/wiki/Betti_number
        https://youtu.be/gVq_xXnwV-4

        """

        raise NotImplementedError()

    @classmethod
    def from_data_frame(cls, df, columns, epsilon, metric=None, prefix=''):
        """

        Create a Vietoris-Rips complex from a pandas DataFrame.

        Parameters:
        -----------
        df: pd.Dataframe
            A dataframe of numeric values.
        columns: list
            List of columns. The columns should contain coordinates of points.
        epsilon: float
            A positive real number.
        metric: callable, optional
            A function that calculates distance between `Point` objects.
            If None the Euclidean metric will be used.
        prefix: str, optional
            The name of a point will be create by concatenating the prefix
            with the index.

        """

        for col in columns:
            if col not in df:
                raise ValueError(f'DataFrame does not contain column {col}')

        pts = list()
        if metric is None:
            metric = euclidean_metric
        for index, row in df[columns].iterrows():
            p = Point(name=f'{prefix}{index}', coords=row.to_list())
            pts.append(p)

        return cls(pts, epsilon, metric)

    @classmethod
    def from_list(cls, names, coords, epsilon, metric=None, prefix=''):
        """

        Create a Vietoris-Rips complex from a list of names and points.

        Parameters:
        -----------
        names: list
            A list of names of points.
        coords: list
            A list of coordinates.
        epsilon: float
            A positive real number.
        metric: callable, optional
            A function that calculates distance between `Point` objects.
            If None the Euclidean metric will be used.
        prefix: str, optional
            A prefix to be added to names.

        """

        if len(names) != len(coords):
            raise ValueError('Labels and coordinates'
                             'should have the same length')

        if metric is None:
            metric = euclidean_metric

        pts = list()
        for name, coord in zip(names, coords):
            p = Point(name=f'{prefix}{name}', coords=coord)
            pts.append(p)
        return cls(pts, epsilon, metric)
