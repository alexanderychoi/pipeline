# -*-coding:Utf-8 -*

import sqlite3
db = sqlite3.connect('data/db_test')
cursor = db.cursor()
cursor.execute('''
    CREATE TABLE users(id TEXT PRIMARY KEY, name TEXT,
                       phone TEXT, email TEXT unique, password TEXT)
''')
db.commit()