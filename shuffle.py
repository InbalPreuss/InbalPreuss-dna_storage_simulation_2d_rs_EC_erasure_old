import sqlite3

from config import PathLike


def shuffle(shuffle_db_file: PathLike, input_file: PathLike, output_file: PathLike):

    shuffle_db_file.unlink(missing_ok=True)
    conn = sqlite3.connect(shuffle_db_file)
    c = conn.cursor()

    c.execute('''CREATE TABLE shuffle_table
                 (idx real, data text)''')

    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            c.execute("INSERT INTO shuffle_table VALUES (" + str(idx) + ", '" + line + "')")

    conn.commit()

    c.execute("SELECT * FROM shuffle_table ORDER BY RANDOM()")
    with open(output_file, 'w+') as f:
        for idx, line in c:
            f.write(line + '\n')


def sort_oligo_file(barcode_len: int, barcode_rs_len: int,
                    sort_db_file: PathLike, input_file: PathLike, output_file: PathLike):
    sort_db_file.unlink(missing_ok=True)
    conn = sqlite3.connect(sort_db_file)
    c = conn.cursor()

    c.execute('''CREATE TABLE sort_table
                     (barcode text, data text)''')

    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            barcode = line[:barcode_len+barcode_rs_len]
            payload = line[barcode_len+barcode_rs_len:]
            c.execute("INSERT INTO sort_table VALUES ('" + barcode + "', '" + payload + "')")

    c.execute("SELECT * FROM sort_table ORDER BY barcode")
    with open(output_file, 'w+') as f:
        for barcode, payload in c:
            f.write(barcode + payload + '\n')
