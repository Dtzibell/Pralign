from PySide6.QtWidgets import QStyleFactory
from typing import Any
from PySide6 import QtCore, QtWidgets, QtGui
import sys
import requests as rq
from ctypes import ArgumentError
import polars as pl
import parasail as ps
import pathlib
import argparse
import re
import os
import tomllib
from .api import get_tm_batch_onetaq_stdbuf, get_tm_batch_phusion_hfbuf


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pralign – primer finder")
        container = QtWidgets.QWidget()
        container.setObjectName("root")
        self.update_button = QtWidgets.QPushButton("Update database")
        self.update_button.clicked.connect(self.magic)
        self.alignment = DBQueryWidget()
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(14)
        layout.addWidget(self.alignment, 1)
        layout.addWidget(self.update_button, 0, QtCore.Qt.AlignmentFlag.AlignRight)
        self.setCentralWidget(container)

    @QtCore.Slot()
    def magic(self):
        answer = QtWidgets.QMessageBox.question(
            self,
            "Update database",
            "Do you really want to update the database? It might take a couple minutes",
        )

        if answer == QtWidgets.QMessageBox.StandardButton.Yes:
            self.update_button.setEnabled(False)
            self._worker = UpdateWorker()
            self._worker.signals.done.connect(self.update_finished)
            self._worker.signals.failed.connect(self.update_failed)
            QtCore.QThreadPool.globalInstance().start(self._worker)

    @QtCore.Slot(object)
    def update_finished(self, _):
        self.update_button.setEnabled(True)
        QtWidgets.QMessageBox.information(self, "Update database", "Database updated.")

    @QtCore.Slot(str)
    def update_failed(self, message):
        self.update_button.setEnabled(True)
        QtWidgets.QMessageBox.critical(self, "Update failed", message)


class UpdateWorker(QtCore.QRunnable):
    def __init__(self):
        super().__init__()
        self.signals = WorkerSignals()

    def run(self):
        try:
            update_db(False)
            self.signals.done.emit(None)
        except Exception as e:
            self.signals.failed.emit(str(e))


class ParamsWidget(QtWidgets.QGroupBox):
    def __init__(self):
        super().__init__("Parameters")
        with open("config.toml", "rb") as f:
            self.defaults = tomllib.load(f)
        d = self.defaults

        self.reverse = QtWidgets.QCheckBox("Reverse complement the query")
        self.reverse.setChecked(d["query"]["reverse"])

        self.top = QtWidgets.QSpinBox()
        self.top.setRange(1, 1000)
        self.top.setValue(d["output"]["primers"]["top"])

        self.shortest = QtWidgets.QSpinBox()
        self.shortest.setMinimum(2)
        self.shortest.setValue(d["filter"]["primers"]["shortest"])
        self.shortest.setSuffix(" nt")

        self.longest = QtWidgets.QSpinBox()
        self.longest.setMinimum(2)
        self.longest.setValue(d["filter"]["primers"]["longest"])
        self.longest.setSuffix(" nt")

        self.output = QtWidgets.QLineEdit(d["output"]["filename"])

        self.shortest.valueChanged.connect(self.pin_longest)
        self.longest.valueChanged.connect(self.pin_shortest)

        form = QtWidgets.QFormLayout(self)
        form.setVerticalSpacing(12)
        form.addRow(self.reverse)
        form.addRow("Top N primers", self.top)
        form.addRow("Shortest", self.shortest)
        form.addRow("Longest", self.longest)
        form.addRow("Output file", self.output)

    @QtCore.Slot()
    def pin_shortest(self):
        if self.shortest.value() > self.longest.value():
            self.shortest.setValue(self.longest.value())

    @QtCore.Slot()
    def pin_longest(self):
        if self.shortest.value() > self.longest.value():
            self.longest.setValue(self.shortest.value())


class WorkerSignals(QtCore.QObject):
    done = QtCore.Signal(object)
    failed = QtCore.Signal(str)


class AlignWorker(QtCore.QRunnable):
    def __init__(self, namespace, defaults):
        super().__init__()
        self.namespace = namespace
        self.defaults = defaults
        self.signals = WorkerSignals()

    def run(self):
        try:
            self.signals.done.emit(run(self.namespace, self.defaults))
        except Exception as e:
            self.signals.failed.emit(str(e))


class DBQueryWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.params = ParamsWidget()

        self.query = QtWidgets.QPlainTextEdit()
        self.query.setPlaceholderText("Enter a DNA sequence")

        self.confirm = QtWidgets.QPushButton("Align primer")
        self.confirm.clicked.connect(self.magic)

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)

        self.table = QtWidgets.QTableWidget()
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectItems
        )
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)

        query_box = QtWidgets.QGroupBox("Query sequence")
        query_box.setMinimumWidth(400)
        q = QtWidgets.QVBoxLayout(query_box)
        q.addWidget(self.query)
        q.addWidget(self.confirm)

        results_box = QtWidgets.QGroupBox("Results")
        results_box.setMinimumWidth(800)
        r = QtWidgets.QVBoxLayout(results_box)
        r.addWidget(self.status)
        r.addWidget(self.table, 1)

        right = QtWidgets.QVBoxLayout()
        right.addWidget(results_box, 1)

        left = QtWidgets.QVBoxLayout()
        left.addWidget(self.params, 0, QtCore.Qt.AlignmentFlag.AlignTop)
        left.addWidget(query_box, 0)

        main = QtWidgets.QHBoxLayout(self)
        main.setContentsMargins(0, 0, 0, 0)
        main.setSpacing(16)
        main.addLayout(left, stretch=1)
        main.addLayout(right, 1)

    @QtCore.Slot()
    def magic(self):
        p = self.params
        namespace = argparse.Namespace(
            query=self.query.toPlainText(),
            reverse=p.reverse.isChecked(),
            top=p.top.value(),
            shortest=p.shortest.value(),
            longest=p.longest.value(),
            output=p.output.text(),
            update=False,
        )
        self.confirm.setEnabled(False)
        self.status.setText("Aligning…")
        self.table.setRowCount(0)

        _worker = AlignWorker(namespace, p.defaults)
        _worker.signals.done.connect(self.show_results)
        _worker.signals.failed.connect(self.show_error)
        QtCore.QThreadPool.globalInstance().start(_worker)

    @QtCore.Slot(object)
    def show_results(self, df: pl.DataFrame):
        self.confirm.setEnabled(True)
        self.status.setText(
            f"{df.height} primers found, saved to {self.params.output.text().strip()}.csv"
        )
        self.table.setColumnCount(len(df.columns))
        self.table.setHorizontalHeaderLabels(df.columns)
        self.table.setRowCount(df.height)
        for i, row in enumerate(df.iter_rows()):
            for j, value in enumerate(row):
                text = f"{value:.3f}" if isinstance(value, float) else str(value)
                item = QtWidgets.QTableWidgetItem(text)
                if not isinstance(value, str):
                    item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()

    @QtCore.Slot(str)
    def show_error(self, message: str):
        self.confirm.setEnabled(True)
        self.status.setText(f"Error: {message}")


def apply_theme(app: QtWidgets.QApplication) -> None:
    app.setStyle(QStyleFactory.create("Fusion"))

    C = QtGui.QColor
    R = QtGui.QPalette.ColorRole
    G = QtGui.QPalette.ColorGroup
    pal = QtGui.QPalette()
    pal.setColor(R.Window, C("#f4f6f9"))
    pal.setColor(R.WindowText, C("#1f2937"))
    pal.setColor(R.Base, C("#ffffff"))  # line edits, spin boxes, table
    pal.setColor(R.AlternateBase, C("#f3f6fa"))  # alternating table rows
    pal.setColor(R.Text, C("#1f2937"))
    pal.setColor(R.PlaceholderText, C("#9aa3ad"))
    pal.setColor(R.Button, C("#ffffff"))
    pal.setColor(R.ButtonText, C("#1f2937"))
    pal.setColor(R.ToolTipBase, C("#1f2937"))
    pal.setColor(R.ToolTipText, C("#ffffff"))
    pal.setColor(R.Highlight, C("#2563eb"))  # selection and focus ring
    pal.setColor(R.HighlightedText, C("#ffffff"))
    pal.setColor(R.Link, C("#2563eb"))
    for role in (R.WindowText, R.Text, R.ButtonText):
        pal.setColor(G.Disabled, role, C("#9aa3ad"))
    app.setPalette(pal)
    app.setStyleSheet(STYLE)


STYLE = """
/* Cards: the title sits on the frame line, as the native style does it */
QGroupBox {
    background: #ffffff;
    border: 1px solid #dfe3e8;
    border-radius: 10px;
    margin-top: 12px;
    padding: 16px 12px 12px 12px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    padding: 0 6px;
    background: #f4f6f9;          /* same as the window, hides the frame line behind the title */
    color: #374151;
    font-weight: 600;
}

/* Only the primary button is custom; all other buttons stay native */
QPushButton#primary {
    background: #2563eb;
    color: white;
    border: none;
    border-radius: 6px;
    padding: 8px 18px;
    font-weight: 600;
}
QPushButton#primary:hover    { background: #1d4ed8; }
QPushButton#primary:pressed  { background: #1e40af; }
QPushButton#primary:disabled { background: #a9c1f5; color: #eef3ff; }

/* Table header only; the table body stays native */
QHeaderView::section {
    background: #f1f3f6;
    color: #374151;
    border: none;
    border-bottom: 1px solid #dfe3e8;
    padding: 6px 8px;
    font-weight: 600;
}
"""


def main() -> None:
    if len(sys.argv) > 1:
        print(f"Running pralign with args: {sys.argv}")
        with open("config.toml", "rb") as f:
            toml = tomllib.load(f)
        argp = argparse.ArgumentParser(
            prog="Pralign",
            description="Find a primer sequence that matches an arbitrarily long query DNA sequence",
        )
        argp.add_argument("query")
        argp.add_argument(
            "--update",
            "-u",
            help="Whether to update the existing database",
            action="store_true",
            default=toml["database"]["update"],
        )
        argp.add_argument(
            "--reverse",
            "-r",
            help="Whether to reverse and make a complement of the query strand",
            action="store_true",
            default=toml["query"]["reverse"],
        )
        argp.add_argument(
            "--top",
            "-n",
            help="The top N primers to output into a csv file",
            default=toml["output"]["primers"]["top"],
        )
        argp.add_argument(
            "--longest",
            "-l",
            help="The length of the longest primer of interest",
            default=toml["filter"]["primers"]["longest"],
        )
        argp.add_argument(
            "--shortest",
            "-s",
            help="The length of the shortest primer of interest",
            default=toml["filter"]["primers"]["shortest"],
        )
        argp.add_argument(
            "--output",
            "-o",
            help='The name of the output file, without the extension. E.g. "alignments"',
            default=toml["output"]["filename"],
        )
        namespace = argp.parse_args()
        run(namespace, toml)
        sys.exit(1)
    else:
        app = QtWidgets.QApplication([])
        apply_theme(app)

        main_window = MainWindow()
        main_window.show()

    sys.exit(app.exec())


# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")
VALID_QUERY = re.compile("^[ACGT]+$")
LABGURU = "https://mftuebingen.labguru.com/api/v1/primers"


def reverse_and_complement(query: str) -> str:
    """Makes a complement and reverses the given query sequence.
    :param query: A query DNA sequence.
    """
    return query.translate(COMP)[::-1]


# TODO: would it not make sense to make gap_open infinity? Would avoid filtering down the line
def align(query: str, seq: str, gap_open: int = 1) -> dict[str, int | str]:
    """Performs Smith-Waterman on the given query and sequence.
    :param query: DNA sequence 1.
    :param seq: DNA sequence 2.
    :param gap_open: penalty for opening and extending a gap.
    """
    result = ps.sw_trace(query, seq, gap_open, gap_open, MATRIX)
    aligned_primer = result.traceback.ref
    return {"alignment score": result.score, "aligned portion": aligned_primer}


def load_dotenv(path=".env"):
    """Minimal .env loader: KEY=value lines, '#' comments, optional quotes.
    Real environment variables take precedence over the file."""
    if not os.path.isfile(path):
        return
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.removeprefix("export ").strip()
            value = value.strip().strip('"').strip("'")
            os.environ.setdefault(key, value)


def update_db(silent_errors: bool) -> None:
    print("Updating the database")
    load_dotenv()
    token = os.environ.get("API_KEY")
    if token is None:
        raise ArgumentError(
            "Could not find API_KEY in .env, hence updating the database is not possible"
        )
    response = True
    counter = 1
    df = pl.DataFrame(
        {"primer name": [], "vial number": [], "sequence": []},
        schema={
            "primer name": pl.String,
            "vial number": pl.String,
            "sequence": pl.String,
        },
    )
    while response:
        print(f"Retrieving page {counter}")
        response = rq.get(
            f"{LABGURU}",
            params={
                "token": token,
                "content_type": "application/json",
                "page": counter,
            },
        ).json()
        counter += 1
        titles = []
        sequences = []
        vial_nos = []
        for primer in response:
            try:
                sequences.append(primer["sequences"][0]["seq"])
                titles.append(primer["title"])
                vial_nos.append(primer["Vial Number"])
            except IndexError as e:
                if not silent_errors:
                    if hasattr(e, "message"):
                        print(
                            f"{e.message}. Likely the primer {primer['title']} contains no sequence"
                        )
                    else:
                        print(
                            f"{e}. Likely the primer {primer['title']} contains no sequence"
                        )
        titles = pl.Series("primer name", titles)
        vial_nos = pl.Series("vial number", vial_nos)
        sequences = pl.Series("sequence", sequences)
        df = pl.concat([df, pl.DataFrame([titles, vial_nos, sequences])])
    df.write_parquet("database.parquet", compression="zstd", compression_level=22)
    print(f"Updated the database. Total primers: {df.height}")


def parse_query(query: str, reverse: bool) -> str:
    query = query.strip().upper()
    if VALID_QUERY.fullmatch(query) is None:
        raise ArgumentError(
            f"query {query} is not valid. Only characters A, C, T, and G are allowed"
        )
    if reverse:
        query = reverse_and_complement(query)
    return query


def parse_shortest(shortest: str) -> int:
    shortest: int = int(shortest)
    if shortest < 2:
        raise ArgumentError(
            f'The argument "shortest" must be at least 2. Given {shortest}.'
        )
    return shortest


def run(namespace: argparse.Namespace, toml: dict[str, Any]) -> pl.DataFrame:
    if namespace.update:
        update_db(toml["runtime"]["silent_errors"])

    query = parse_query(namespace.query, namespace.reverse)
    top: int = int(namespace.top)
    longest: int = int(namespace.longest)
    shortest = parse_shortest(namespace.shortest)
    output: str = namespace.output.strip()

    print("Aligning primers")
    db: pathlib.Path = pathlib.Path("database.parquet")
    df: pl.DataFrame = (
        pl.scan_parquet(db)
        .with_columns(pl.col("sequence").str.len_bytes().alias("sequence length"))
        .filter(
            (pl.col("sequence length") < longest)
            & (pl.col("sequence length") > shortest)
        )  # some sequences are question marks? Not relevant, but must set lower bound to at least 2
        .with_columns(
            pl.col("sequence")
            .map_elements(
                lambda seq: align(query, seq.upper()),
                return_dtype=pl.Struct(
                    {"alignment score": pl.Int64, "aligned portion": pl.String}
                ),
            )
            .alias("alignment")
        )  # performs the alignment
        .unnest("alignment")  # unnests the struct into two separate columns
        .filter(
            pl.col("aligned portion")
            .str.contains("-")
            .not_()
            .and_(pl.col("alignment score") > 15)
        )  # filters out gapped and insufficient alignments. TODO: introduce a toml for params? Overkill?
        .sort(by="alignment score", descending=True)
        .with_columns(
            (pl.col("alignment score") / pl.col("sequence length")).alias(
                "aligned proportion"
            )
        )  # Computes the aligned proportion of the primer. Alignment score = alignment length because match score is 1.
        .sort(by="aligned proportion", descending=True)
        .head(top)
        .collect()
    )
    onetaq_tm = pl.Series(
        "OneTaq Tm",
        [
            dct["tm1"]
            for dct in get_tm_batch_onetaq_stdbuf(df["aligned portion"].to_list())
        ],
    )
    phusion_tm = pl.Series(
        "Phusion Tm",
        [
            dct["tm1"]
            for dct in get_tm_batch_phusion_hfbuf(df["aligned portion"].to_list())
        ],
    )
    df = df.with_columns(onetaq_tm, phusion_tm).drop("aligned portion")
    df.write_csv(f"{output}.csv", float_precision=3)
    print(f"Alignment output written to {output}.csv")
    return df
