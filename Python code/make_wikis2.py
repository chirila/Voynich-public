"""
August 13, 2019
Runs make_wiki_corpus_short.py on all files in a subdirectory
"""

import sys
import re
from os import listdir
from gensim.corpora import WikiCorpus

def make_corpus(in_f, out_f, num_f):
    """Convert Wikipedia xml dump file to text corpus"""

    output = open(out_f, 'w')
    wiki = WikiCorpus(in_f, token_min_len=1, token_max_len=5000, lower=False)

    print('Corpus Created')
    full_text = list(wiki.get_texts())
    i = 0
    for text in full_text[0:int(num_f)]:

        output.write(bytes(' '.join(text).encode('utf-8')) + '\n')
        i = i + 1
        if (i % 100 == 0):
            print('Processed ' + str(i) + ' articles')

    output.close()
    print('Processing complete!')


def run_wikicorpus(in_dir, num_articles):

    for f in listdir(in_dir):
        
        if str(f).endswith('bz2'):
            print(f)

            a =  in_dir + '/' + f
            b = in_dir + '/' + re.sub(r'\-.*$', '.txt', f)

            make_corpus(a, b, num_articles)

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('Usage: python test.py <directory> <article length>')
        sys.exit(1)
    in_dir = sys.argv[1]
    num_articles = sys.argv[2]

    run_wikicorpus(in_dir, num_articles)