"""
SeqLit bioinformatics modules package.
Configures system-wide SSL certificates using certifi to ensure
reliable HTTPS connectivity to NCBI E-utilities, NCBI BLAST, and UniProt REST APIs.
"""
import ssl

try:
    import certifi
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except ImportError:
    pass
