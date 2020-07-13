import sqlite3
from sqlite3 import Error

def create_connection(db_file: str) -> sqlite3.Connection:
    """ create a database connection to a SQLite database """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn

def create_table(conn: sqlite3.Connection, create_table_sql: str) -> bool:
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
        return True

    except Error as e:
        print(e)
        return False

def execute_indel(conn: sqlite3.Connection, statement: str) -> bool:
    '''
    Execute a command. Returns True if the command could be executed, False otherwise

    Inputs:
        conn:       (sqlite3.Connection) the connection to the database
        statement:  (str) the sql statement to use for insert or deletion
    Outputs:
        (bool) True if successfully executed, False otherwise
    '''
    try:
        c = conn.cursor()
        c.execute(statement)
        conn.commit()
        return True

    except Error as e:
        print(f'ERROR in execute_indel:\n{e}')
        print(f'SQL query was:\n{statement}')
        return False

def execute_indel_many(conn: sqlite3.Connection, statement: str, many: list) -> bool:
    '''
    Execute many commands. Returns True if the commands could be executed, False otherwise

    Inputs:
        conn:       (sqlite3.Connection) the connection to the database
        statement:  (str) the statement to run. It should have ? in places where the tuples of many go 
                            (see https://docs.python.org/3/library/sqlite3.html#sqlite3.Connection.executemany)
        many:       (list) the many commands to excute
    Outputs:
        (bool) True if successfully executed, False otherwise
    '''
    try:
        c = conn.cursor()
        c.executemany(statement, many)
        conn.commit()
        return True

    except Error as e:
        print(f'ERROR in execute_indel_mane:\n{e}')
        print(f'SQL query was:\n{statement}')
        return False

def execute_return(conn: sqlite3.Connection, statement: str) -> sqlite3.Cursor:
    '''
    Execute a command that has a return

    Inputs:
        conn:       (sqlite3.Connection) the connection to the database
        statement:  (str) the sqlite statement to execute
    Outputs:
        (sqlite3.Cursor) the return value of the query
    '''
    try:
        c = conn.cursor()
        return c.execute(statement)

    except Error as e:
        print(f'ERROR in execute_return:\n{e}')
        print(f'SQL query was:\n{statement}')
        return None

