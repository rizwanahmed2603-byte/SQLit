"""
SQLite caching module for SeqLit.
Prevents redundant BLAST, Entrez, and UniProt requests by caching results by sequence MD5 hash.
"""

import sqlite3
import hashlib
import json
import os
from typing import Optional, Dict, Any

DB_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "seqlit_cache.db")

def get_connection():
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS analysis_cache (
            sequence_hash TEXT PRIMARY KEY,
            sequence_type TEXT,
            sequence_length INTEGER,
            result_json TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    conn.commit()
    conn.close()

def compute_hash(sequence: str) -> str:
    return hashlib.md5(sequence.strip().upper().encode("utf-8")).hexdigest()

def get_cached_result(sequence: str) -> Optional[Dict[str, Any]]:
    seq_hash = compute_hash(sequence)
    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT result_json FROM analysis_cache WHERE sequence_hash = ?", (seq_hash,))
    row = cursor.fetchone()
    conn.close()
    if row:
        try:
            return json.loads(row["result_json"])
        except Exception:
            return None
    return None

def set_cached_result(sequence: str, seq_type: str, result: Dict[str, Any]):
    seq_hash = compute_hash(sequence)
    conn = get_connection()
    cursor = conn.cursor()
    cursor.execute("""
        INSERT OR REPLACE INTO analysis_cache (sequence_hash, sequence_type, sequence_length, result_json)
        VALUES (?, ?, ?, ?)
    """, (seq_hash, seq_type, len(sequence), json.dumps(result)))
    conn.commit()
    conn.close()

# Auto-initialize database on import
init_db()
