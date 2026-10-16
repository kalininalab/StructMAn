"""
Central MySQL driver selection for StructMan.

Bulk runs issue very large multi-row INSERT statements (tens of millions of rows
per call for whole-PDB inputs). In that path the dominant client-side cost is
escaping every parameter, and pymysql does that in pure Python. mysqlclient is a
thin wrapper around libmysqlclient and does the same work in C, so we prefer it
whenever it is importable and fall back to pymysql otherwise.

StructMan only uses the portable DB-API 2.0 surface (connect, cursor, execute,
fetchone, fetchall, commit, rollback, close), so the two drivers are
interchangeable apart from the connect-time differences handled below.

Install the faster driver with:  pip install mysqlclient
(needs libmysqlclient-dev / mysql-devel and a compiler at install time)
"""

import os

DRIVER_NAME: str

# STRUCTMAN_MYSQL_DRIVER forces a driver instead of preferring mysqlclient.
# Useful when mysqlclient is installed but unusable in the current environment
# (a conda build missing its auth plugins, say), and for A/B timing the two
# without touching the environment.
_forced = os.environ.get('STRUCTMAN_MYSQL_DRIVER', '').strip().lower()

if _forced == 'pymysql':
    import pymysql as _driver
    DRIVER_NAME = 'pymysql'
elif _forced == 'mysqlclient':
    import MySQLdb as _driver
    DRIVER_NAME = 'mysqlclient'
elif _forced:
    raise ValueError(
        f"STRUCTMAN_MYSQL_DRIVER={_forced!r} is not a known driver; "
        "use 'mysqlclient' or 'pymysql'")
else:
    try:
        import MySQLdb as _driver
        DRIVER_NAME = 'mysqlclient'
    except ImportError:
        import pymysql as _driver
        DRIVER_NAME = 'pymysql'


# DB-API exceptions, re-exported so callers never need to know the driver.
Error = _driver.Error
Warning = _driver.Warning
InterfaceError = _driver.InterfaceError
DatabaseError = _driver.DatabaseError
DataError = _driver.DataError
OperationalError = _driver.OperationalError
IntegrityError = _driver.IntegrityError
InternalError = _driver.InternalError
ProgrammingError = _driver.ProgrammingError
NotSupportedError = _driver.NotSupportedError

# pymysql exposes these under `.err` as well; keep that spelling working so
# `except mysql_driver.err.OperationalError` stays valid either way.
err = getattr(_driver, 'err', _driver)


def connect(host=None, user=None, password=None, database=None, autocommit=False, **kwargs):
    """
    Driver-agnostic connect.

    Positional order matches the historical `MySQLdb.connect(host, user, pw, db)`
    call sites. Note that pymysql >= 1.0 made all of its connect arguments
    keyword-only, so those positional calls only work through this wrapper.
    """
    # Only options both drivers accept are passed here. mysqlclient rejects
    # unknown keywords outright (max_allowed_packet, for instance, is pymysql
    # only), and its client-side packet limit is not settable through connect()
    # anyway -- it defaults to 1GB in libmysqlclient, well above the server's
    # max_allowed_packet that Config.max_package_size already sizes against.
    params = {
        # pymysql >= 1.0 defaults to utf8mb4, while mysqlclient inherits the
        # charset from the server/my.cnf. Pin it so both drivers decode text
        # columns identically; that also makes both return str for text columns
        # without needing use_unicode, which is the default on both. Binary and
        # BLOB columns still come back as bytes, which the packed-alignment code
        # relies on.
        'charset': 'utf8mb4',
        'autocommit': autocommit,
    }

    if host is not None:
        params['host'] = host
    if user is not None:
        params['user'] = user
    if password is not None:
        params['password'] = password
    if database is not None:
        params['database'] = database

    params.update(kwargs)

    return _driver.connect(**params)
