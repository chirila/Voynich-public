"""
Creates a corpus from Wikipedia dump file.
Inspired by:
https://github.com/panyang/Wikipedia_Word2vec/blob/master/v1/process_wiki.py

Luke Lindemann found this code snippet at 11/29/18 on
https://www.kdnuggets.com/2017/11/building-wikipedia-text-corpus-nlp.html

Must be specified at command line, i.e. python make_wiki_corpus.py enwiki-latest-pages-articles.xml.bz2 wiki_en.txt


### CHANGED: Breaks after num_f many files processed
### CHANGED Default minimum token is 2, maximum is 15, changed to 1,30
### CHANGED No default lowercase
"""

import sys
from gensim.corpora import WikiCorpus

def make_corpus(in_f, out_f, num_f):
    """Convert Wikipedia xml dump file to text corpus"""

    output = open(out_f, 'w')
    wiki = WikiCorpus(in_f, token_min_len=1, token_max_len=50, lower=False)

    print ('Corpus Created')

    i = 0
    for text in wiki.get_texts():
        output.write(bytes(' '.join(text).encode('utf-8')) + '\n')
        i = i + 1
        if (i % 100 == 0):
            print('Processed ' + str(i) + ' articles')
        if (i >= int(num_f)):
            break
    print('Closing...')
    output.close()
    print('Processing complete!')



if __name__ == '__main__':

    if len(sys.argv) != 4:
        print('Usage: python make_wiki_corpus.py <wikipedia_dump_file> <processed_text_file> <num_articles_to_process>')
        sys.exit(1)
    in_f = sys.argv[1]
    out_f = sys.argv[2]
    num_f = sys.argv[3]
    make_corpus(in_f, out_f, num_f)