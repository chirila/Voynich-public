#### This is a suite of functions for Entropy Statistics
#### It uses the 'ngram' and 'stringdist' packages for ngram calculation
#### 'lindentropy' takes a string and summarizes entropy character statistics for 1-3 orders
#### 'condbigram' and 'condtrigram' take a string and give a summary of bigram/trigram absolute and conditional
####      probabilities
### Luke Lindemann April 17, 2019

rm(list=ls())
library(ngram)
library(stringdist) # ngram doesn't work with non-Latin characters, use 'qgrams'
library(stringr)


# Here is a function that summarizes the entropy profile of a character string. 
#   Entropy is expressed in bits (shannons) 
#   h0 = Zero Order Entropy = maximum entropy given the num of characters
#   H1 = First Order Absolute Entropy (also h1) = character entropy
#   H2 = Second Order Absolute Entropy = digraph entropy
#   h2 = Second Order Conditional Entropy = H2-H1 =


sumentropy <- function (s, remove.spaces = FALSE)
{
  if (nchar(s) <= 1)
  {
    return(NA)
  } 
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }

  s.list = strsplit(s, "")
  s.freq = table(s.list)
  char.num = length(s.freq)
  s.norm <- s.freq/sum(s.freq)
  
  s.freq.bi = qgrams(s, q=2)
  s.norm.bi <- s.freq.bi/sum(s.freq.bi)
  
  s.freq.tri = qgrams(s, q=3)
  s.norm.tri <- s.freq.tri/sum(s.freq.tri)
  
  h0 = log2(length(s.freq))
  
  H1 = -sum(log2(s.norm)*s.norm)
  
  H2 = -sum(log2(s.norm.bi)*s.norm.bi)
  
  H3 = -sum(log2(s.norm.tri)*s.norm.tri)
  
  h2 = H2-H1
  
  h3 = H3-H2
  
  diff = H1-h2 
  
  quot = h2/h0
  
  s.summary = as.table(c(char.num, h0,H1,H2,h2,diff,quot, H3, h3))
  rownames(s.summary) = c("char.num","h0", "h1", "H2", "h2", "h1-h2", "h2/h0", "H3", "h3")
  s.summary
  
}



# Creates a chart of Bigrams from a Character String:
  
condbigram <- function (s, remove.spaces = FALSE) {
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }
  
  char.freq <- (as.data.frame(t(qgrams(freq=s, q=1))))
  char.freq$ngrams <- rownames(char.freq)
  char.freq$prop <- char.freq$freq / sum(char.freq$freq)
  bi.freq <- (as.data.frame(t(qgrams(freq=s, q=2))))
  bi.freq$ngrams <- rownames(bi.freq)
  bi.freq$prop <- bi.freq$freq / sum(bi.freq$freq)
  
  cond.Prop <- c()
  fchar.Prop <- c()
  
  
  for (gram in bi.freq$ngrams) {
    
    # First Character Freq
    charfreq <- char.freq[char.freq$ngrams == substr(gram,1,1), ]$freq
    
    # First Character Prop
    charprop <- char.freq[char.freq$ngrams == substr(gram, 1, 1), ]$freq / sum(char.freq$freq)
    
    # Conditional Freq = Bigram Freq / First Char Frequency
    freqdiv <- bi.freq[bi.freq$ngrams == gram,]$freq / charfreq
    
    fchar.Prop <- c(fchar.Prop, charprop)
    cond.Prop <- c(cond.Prop, freqdiv)
    
  }
  
  bi.freq$cond.prop <- cond.Prop 
  bi.freq$fchar.prop <- fchar.Prop
  bi.freq[order(bi.freq$cond.prop, decreasing = TRUE),]
}




# Predictive Conditional Bigrams: a chart of all bigrams with >50% conditional probability, ranked by absolute probability

predcondbigram <- function (s, rm.spaces = FALSE) {
  
  df <- condbigram(s, remove.spaces = rm.spaces)  
  
  df <- df[df$cond.prop >= 0.5,]
  
  df[order(df$prop, decreasing=TRUE),]
  
}




condtrigram <- function (s, remove.spaces = FALSE) {
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }
  
  bi.freq <- (as.data.frame(t(qgrams(freq=s, q=2))))
  bi.freq$ngrams <- rownames(bi.freq)
  bi.freq$prop <- bi.freq$freq / sum(bi.freq$freq)
  tri.freq <- (as.data.frame(t(qgrams(freq=s, q=3))))
  tri.freq$ngrams <- rownames(tri.freq)
  tri.freq$prop <- tri.freq$freq / sum(tri.freq$freq)
  
  
  cond.Prop <- c()
  pair.freq <- c()
  
  
  for (gram in tri.freq$ngrams) {
    
    # First two characters
    pair <- substring(gram, 1,2)
    # Frequency of first two characters
    pair.freq <- bi.freq[bi.freq$ngrams == pair,]$freq
    # Conditional probability: frequency of trigram / frequency of character pair
    freqdiv <- tri.freq[tri.freq$ngrams == gram,]$freq / pair.freq
    
    cond.Prop <- c(cond.Prop, freqdiv)
    
  }
  
  tri.freq$cond.prop <- cond.Prop
  
  tri.freq[order(tri.freq$cond.prop, decreasing = TRUE),]
}



# This is a homemade bigram function that is less efficient than ngram and stringdist

#lindebigram <- function (s, remove.spaces = FALSE)
# WORKS FOR SHORT MESSAGES
#{
# Take out spaces if remove.spaces = TRUE
#if (remove.spaces == TRUE)
#{
# s = gsub(" ", "", s)
#} 

#s.list <- strsplit(s, "")

# A data frame of individual character frequency
#char.freq <- as.data.frame(table(s.list))


# A list of all bigrams
#bi.list <- c()

#for (char in char.freq$s.list) {

# char.list <- regmatches (s, gregexpr(paste (char, ".", sep = ""), s, perl = TRUE ) )

#bi.list <- c(bi.list, unlist(char.list))

#}

# A frequency table and prop table of all bigrams
#bigrams <- as.data.frame (table(bi.list))

#bigrams$Prop <- prop.table(bigrams$Freq)

# A conditional frequency table for all bigrams
#cond.Prop <- c()

#for (gram in bigrams$bi.list) {

# The character frequency for the character that begins the bigram
# charfreq <- char.freq[char.freq$s.list == substr(gram, 1,1), ]$Freq

# The bigram frequency divided by the character frequency
#freqdiv <- bigrams[bigrams$bi.list == gram,]$Freq / charfreq

#cond.Prop <- c(cond.Prop, freqdiv)
#}

#bigrams$cond.Prop <- cond.Prop

#bigrams[order(bigrams$Prop, decreasing = TRUE),]

#}




### THIS FUNCTION GOES THROUGH EACH STEP OF THE SUKHOTIN ALGORITHM FOR DETERMINING
### THE VOWELS IN AN ALPHABETIC CHARACTER STRING

sukhotin_vowels <- function (s, remove.spaces = FALSE, show.steps = FALSE) {
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }

  # Initialize matrix
  s.bi <- qgrams(s, q=2)
  s.char <- table(unlist(strsplit(s,"")))
  s.mat <- matrix (0L, nrow=length(s.char), ncol=length(s.char)+1)
  rownames(s.mat) <- names(s.char)
  colnames(s.mat) <- c(names(s.char), "Sum")
  s.cv <- rep('c', length(s.char))
  names(s.cv) <- names(s.char)
  
  # Populate with adjacent values using the bigrams (excluding the diagonal)
  # (Steps 1-2)
  for (r in names(s.char)) {
    for (c in names(s.char)) {
      
      if (!(r == c)) {
        rc = paste(r,c,sep="")
        cr = paste(c,r,sep="")
        adj_ct = 0
        
        if (rc %in% colnames(s.bi)) {
          adj_ct = adj_ct + s.bi[,rc]
        }
        if (cr %in% colnames(s.bi)) {
          adj_ct = adj_ct + s.bi[,cr]
        }       
        
        s.mat[r,c] = adj_ct
      }
    }
  }
  
  # Add sums (Step 3)
  for (r in names(s.char)) {
    s.mat[r,'Sum'] = sum(s.mat[r,])
  }

  if (show.steps) {
    print ("Initial Matrix:")
    print(s.mat)
    print(s.cv)
  }
  
  foundVowels = FALSE
  
  while (foundVowels == FALSE) {
   
    # Find the consonant with the highest sum. Call that a vowel.
    # (Step 4)
    
    # This is a version of the matrix where all vowels are given a sum of 0 so that only consonants are 
    # compared
    s.mat.cons <- s.mat
    for (chr in names(s.char)) {
      if (s.cv[chr] == 'v') {
        s.mat.cons[chr,'Sum'] = 0
      }
    }
    
    
    high_char <- names(which.max(s.mat.cons[,'Sum']))
    
    # If the highest value is still greater than 0, don't terminate
    if (s.mat.cons[high_char,'Sum'] > 0) {

      # This is the new vowel
      s.cv[high_char] <- 'v'
      
      if (show.steps) {
        print("New Vowel: ")
        print(high_char)
        print(s.cv)
      }

      
      
      # For each consonant, find its adjacency count to the new vowel.
      # Subtract twice that amount from the consonant's sum
      # (Step 5)
      
      for (chr in names(s.char)) {
        
        if (s.cv[chr] == 'c') {
          adj_ct <- s.mat[chr,high_char]
          s.mat[chr,'Sum'] = s.mat[chr,'Sum'] - (2 * adj_ct)
        }
      }
      
      if (show.steps) {
        print('Adjusted Matrix')
        print(s.mat)
        #readline(prompt="Press [enter] to continue")
      }
      
    } else {
      
    foundVowels = TRUE
    print('Vowels: ')
    for (chr in names(s.char)) {
      
      if (s.cv[chr] == 'v' ) {
        
        print (chr)
      }
    }
   }
  }
}





## Given an index on a word list, pulls out consecutive repetitions

rep.list <- function (x, wlist, distance, wordmin) {
  
  if  (x < (length(wlist)))   {
    
    if (nchar(wlist[x]) >= wordmin) {
      if (adist(wlist[x],wlist[x+1]) <= distance) {
        return(c(wlist[x], rep.list(x+1, wlist, distance, wordmin)  )  )
        
      } else {
        
        return (c(wlist[x]))  
      }
    }
  } else {
    
    return (c(wlist[x]))  
    
  }
  
}




## This function pulls out consecutive repetitive sequences of a text
  # dist is the levenshtein distance
  # min.seq is the minimum length of a sequence
  # min.word is the minimum length of a word

rep.seqs <- function (s, dist = 2, min.seq = 2, min.word = 3) {
  
  wlist <- unlist(strsplit(s, " "))
  
  i = 1
  while (i <= length(wlist)) {
    
    seq <- rep.list(i, wlist, dist, min.word)

    if (length(seq) > 1) {
      
      if (length(seq) >= min.seq) {
        
        print(seq)        
      
      }
      i = i + length(seq) - 1
      
    } else {
     
      i = i + 1 
    }
    
  }
  
}



### This function computes the MATTR (Moving Average Type Token Ratio)
  # for a document. The Type-Token ratio is computed at different windows
  # and then averaged

mattr <- function(s, window = 500) {
  
  word.list <- unlist(strsplit(s, " "))
  
  # Super-useful partition function for R:
  word.part <- split(word.list, as.integer((seq_along(word.list) - 1) / window))
  
  # Type-token frequency for each partion
  ttr.list <- c()
  for (x in word.part) {
    
    token.freq <- length(x)
    type.freq <- length(table(x))
    ttr.list <- c(ttr.list, type.freq/token.freq)
  }
  
  # Return the avg
  return (mean(ttr.list))
  
}


### This function produces a table of the most common initial or final segments in a text
# Useful for looking at prefixes and suffixes
#   final = FALSE (prefixes), final = TRUE (suffixes)
#   segs = number of segments

affixes <- function(s, final = FALSE, segs = 2) {

  word.list <- unlist(strsplit(s, " "))

  if (final) {
    segments <- str_sub(word.list, -segs, -1)

  } else {
    segments <- str_sub(word.list, 1, segs)

  }

  seg.freq <- as.data.frame(table(segments))

  seg.freq$Prop <- seg.freq$Freq / sum(seg.freq$Freq)

  seg.freq[order(seg.freq$Freq, decreasing = TRUE),]

}


### This function takes a string and returns the string reversed

rev_str <- function(s) {
  return (paste(rev(unlist(strsplit(s,""))),collapse=""))
  
}

  

#### This compares the differences in relative character frequencies between two texts. 

compare_char_freqs <- function(base_text, comp_text, absolute = TRUE) {
  base.text.chars <- unlist(strsplit(base_text, ""))
  comp.text.chars <- unlist(strsplit(comp_text, ""))
  
  full.char.list <- union(base.text.chars, comp.text.chars)
  
  base.text.freqs <- table(factor(base.text.chars, levels = full.char.list))
  comp.text.freqs <- table(factor(comp.text.chars, levels = full.char.list))
  
  base.text.rel.freqs <- base.text.freqs / (sum(base.text.freqs))
  comp.text.rel.freqs <- comp.text.freqs / (sum(comp.text.freqs))
  
  if (absolute) {
    multiplier = 1/100
  } else {
    multiplier = comp.text.rel.freqs
  }
  
  
  comp.table <- (comp.text.rel.freqs - base.text.rel.freqs) / multiplier
  
  comp.table[order(comp.table, decreasing = TRUE)]
  
  
}




#### This compares the differences in relative word frequencies between two texts. 

compare_freqs <- function(base_text, comp_text, sep="") {
  base.text.words <- unlist(strsplit(base_text, sep))
  comp.text.words <- unlist(strsplit(comp_text, sep))
  
  full.word.list <- union(base.text.words, comp.text.words)
  
  base.text.freqs <- table(factor(base.text.words, levels = full.word.list))
  comp.text.freqs <- table(factor(comp.text.words, levels = full.word.list))
  
  base.text.rel.freqs <- base.text.freqs / (sum(base.text.freqs))
  comp.text.rel.freqs <- comp.text.freqs / (sum(comp.text.freqs))

  
  comp.table <- (comp.text.rel.freqs / base.text.rel.freqs)
  
  comp.table[order(comp.table, decreasing = TRUE)]
  
  
}



#### This compares the differences in relative n-gram frequencies between two texts. 

compare_grams <- function(base_text, comp_text, n=2) {
  
  base.text <- gsub(' ', '\\#', base_text)
  comp.text <- gsub(' ', '\\#', comp_text)
  
  base.text.freqs <- qgrams(base.text, q=n)
  comp.text.freqs <- qgrams(comp.text, q=n)
  
  base.text.rel.freqs <- base.text.freqs / sum(base.text.freqs)
  comp.text.rel.freqs <- comp.text.freqs / sum(comp.text.freqs)
  
  base.text.df <- as.data.frame(t(base.text.rel.freqs))
  comp.text.df <- as.data.frame(t(comp.text.rel.freqs))
  
  full.text.df <- merge(base.text.df, comp.text.df, by = 0, all = TRUE)
  
  colnames(full.text.df) <- c('ngram', 'base.freq', 'comp.freq')

  full.text.df$comp <- full.text.df$comp.freq / full.text.df$base.freq
  
  full.text.df[order(full.text.df$comp.freq, decreasing=TRUE),]
  
}








##### Measures of bigram Complexity
# Conditional Character Entropy
# Percentage of bigram space covered
# Standard deviation of bigram percentages
# Perplexity


bicomp <- function (s, remove.spaces = FALSE)
{
  if (nchar(s) <= 1)
  {
    return(NA)
  } 
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }
  
  s.list = strsplit(s, "")
  s.freq = table(s.list)
  char.num = length(s.freq)
  s.norm <- s.freq/sum(s.freq)
  
  s.freq.bi = qgrams(s, q=2)
  s.norm.bi <- s.freq.bi/sum(s.freq.bi)
  
  s.freq.tri = qgrams(s, q=3)
  s.norm.tri <- s.freq.tri/sum(s.freq.tri)
  
  # Conditional Entropy
  H1 = -sum(log2(s.norm)*s.norm)
  H2 = -sum(log2(s.norm.bi)*s.norm.bi)
  h2 = H2-H1
  
  # Predictive Bigram Percentage
  doc.bi <- condbigram(s, remove.spaces=FALSE)
  doc.bi.pred <- doc.bi[doc.bi$cond.prop >= 0.5,] 
  pred.bi.perc <- sum(doc.bi.pred$prop) * 100
 
  # Bigram Space Covered
  s.bi.poss.num <- char.num^2
  bi.cover <- length(s.norm.bi) / s.bi.poss.num * 100
  
  # Standard Deviation of Bigram Proportions
  bi.props <- sd(c(s.norm.bi, rep(0, s.bi.poss.num-length(s.norm.bi))))
  
  # Perplexity
  perp <- 2^(H2) # Should this be H2 or h2?

  s.summary = as.table(c(round(h2,3), round(pred.bi.perc,3), round(bi.cover,3), round(bi.props,3), round(perp,3)))
  rownames(s.summary) = c("Conditional Entropy","Predictive Bigram Percentage", "Bigram Space Covered", "Bigram Std Dev", "Perplexity")
  s.summary
  
}



##### Percentage of Bigrams Coverage for different proportions of predictability
bicovperc <- function (s, perc, remove.spaces = FALSE) {

  doc.bi <- condbigram(s, remove.spaces=remove.spaces)
  
  doc.bi.perc <- doc.bi[doc.bi$cond.prop >= perc,] 
  
  return (sum(doc.bi.perc$prop) / sum(doc.bi$prop))
  
}


### Cosine Similarity Function for two vectors
cos.sim <- function(A, B) 
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}





### Word Entropy
  # Note: switching to ngram bc I think it works for non-Latin at the word level
  # qgram only does character-level stats

wordentropy <- function (s)
{
  if (nchar(s) <= 1)
  {q
    return(NA)
  } 
  
  
  s.list = strsplit(s, " ")
  s.freq = table(s.list)
  word.num = length(s.freq)
  s.norm <- s.freq/sum(s.freq)
  
  s.freq.bi.ngram = ngram(s, n=2)
  s.freq.bi = get.phrasetable(s.freq.bi.ngram)
  s.norm.bi <- s.freq.bi$prop
  
  s.freq.tri.ngram = ngram(s, n=3)
  s.freq.tri = get.phrasetable(s.freq.tri.ngram)
  s.norm.tri <- s.freq.tri$prop
  
  h0 = log2(length(s.freq))
  
  H1 = -sum(log2(s.norm)*s.norm)
  
  H2 = -sum(log2(s.norm.bi)*s.norm.bi)
  
  H3 = -sum(log2(s.norm.tri)*s.norm.tri)
  
  h2 = H2-H1
  
  h3 = H3-H2
  
  diff = H1-h2 
  
  quot = h2/h0
  
  s.summary = as.table(c(word.num, h0,H1,H2,h2,H3,h3))
  rownames(s.summary) = c("word.num","h0", "h1", "H2", "h2", "H3", "h3")
  s.summary
  
}



### Neighbor List : if r=TRUE, the input is a regular expression, otherwise it is a word

nbr.list <- function (s, w, dist = -1, r=TRUE)
{
  
  if (!r) {g
    w <- paste('^', w, '$', sep='')
  }
  
  s.list <- unlist(strsplit(s, " "))
  
  nbrs <- c()
  
  for (i in 1:(length(s.list)))  {
    
    if (grepl(w, s.list[i])) {
      
      if (dist + i <= length(s.list)) {
        nbrs <- c(nbrs, s.list[dist+i])  
      }
    }
  }
  
  return (nbrs)
}



#### Word Frequency Chart --- (I keep having to write this code)


word.freq <- function (s) {
 
  s.words <- unlist(strsplit(s, " "))
  
  s.freq <- as.data.frame(table(s.words))
  
  s.freq$Prop <- s.freq$Freq / sum(s.freq$Freq)
  
  s.freq <- s.freq[order(s.freq$Freq, decreasing = TRUE),]
  
  s.freq$Rank <- c(1:length(s.freq$s.words))
  
  s.freq
}





### Immediate repetitiveness: What percentage of word bigrams are repetitive (dist = Levenshtein distance)

word.rep <- function (s, dist=0) {
  bigrams <- ngram(s, n=2)
  bi.table <- get.phrasetable(bigrams)
  nlist <- c()
  for (n in bi.table$ngrams) {
    if (adist(str_extract(n, '^[^ ]*'), str_sub(str_extract(n, '[^ ]* $'), 1, -2)) <= dist) {
      nlist <- c(nlist, n)
    }
  }
  rep.bi.table <- bi.table[bi.table$ngrams %in% nlist,]
  
  return(sum(rep.bi.table$prop))
}






### Find just the h2 rather than the whole sumentropy summary

quick_h2 <- function (s, remove.spaces = FALSE)
{
  if (nchar(s) <= 1)
  {
    return(NA)
  } 
  
  # How to treat spaces
  if (remove.spaces){
    s <- gsub(' ', '', s)
  }else{
    s <- gsub(' ', '\\#', s)
  }
  
  s.list = strsplit(s, "")
  s.freq = table(s.list)
  s.norm <- s.freq/sum(s.freq)
  
  s.freq.bi = qgrams(s, q=2)
  s.norm.bi <- s.freq.bi/sum(s.freq.bi)
  
  
  H1 = -sum(log2(s.norm)*s.norm)
  
  H2 = -sum(log2(s.norm.bi)*s.norm.bi)
  
  h2 = H2-H1
  
  return(h2)
}





### Delta Entropy: Chart of how much each bigram contributes to conditional character entropy (h2)

delta_h2 <- function (s, remove.spaces=FALSE) {
  
  # Remove spaces or replace spaces with '#'
  if (remove.spaces) {
    s <- gsub(' ', '', s)
  } else {
    s <- gsub(' ', '#', s)
  }
  
  # Asterisks are regexes so I'm replacing them with 'U'
  s <- gsub('\\*', 'U', s)
  
  # h2 of the unaltered text
  base_h2 <- quick_h2(s)
  
  # List of all bigrams and their frequencies
  bigram_table <- qgrams(s, q=2)
  bigram <- colnames(bigram_table)
  freq <- as.vector(bigram_table/sum(bigram_table)) # proportional frequency of each bigram
  #print(length(bigram))
  

  # Change in h2 when each bigram is deleted
  
  i = 0
  bi_delta_h2 <- c()
  for (b in bigram) {
    new_s <- gsub(b, '', s)
    new_h2 <- quick_h2(new_s)
    bi_delta_h2 <- c(bi_delta_h2, new_h2 - base_h2)
    #print(i/length(bigram)*100)
    i = i + 1
  }
  
  # Return sorted dataframe
  df <- data.frame(bigram, freq, bi_delta_h2)
  return(df[order(df$bi_delta_h2, decreasing=TRUE),])
  
}






delta_h3 <- function (s, remove.spaces=FALSE) {
  
  # Remove spaces or replace spaces with '#'
  if (remove.spaces) {
    s <- gsub(' ', '', s)
  } else {
    s <- gsub(' ', '#', s)
  }
  
  # Asterisks are regexes so I'm replacing them with 'U'
  s <- gsub('\\*', 'U', s)
  
  # h3 of the unaltered text
  base_h3 <- sumentropy(s)[8]
  
  # List of all trigrams and their frequencies
  trigram_table <- qgrams(s, q=3)
  trigram <- colnames(trigram_table)
  freq <- as.vector(trigram_table/sum(trigram_table)) # proportional frequency of each trigram
  #print(length(trigram))
  
  
  # Change in h3 when each trigram is deleted
  
  i = 0
  tri_delta_h3 <- c()
  for (t in trigram) {
    new_s <- gsub(t, '', s)
    new_h3 <- sumentropy(new_s)[8]
    tri_delta_h3 <- c(tri_delta_h3, new_h3 - base_h3)
    print(i/length(trigram)*100)
    i = i + 1
  }
  
  # Return sorted dataframe
  df <- data.frame(trigram, freq, tri_delta_h3)
  return(df[order(df$tri_delta_h3, decreasing=TRUE),])
  
}





zipf_fit <- function (s, end.cutoff = 100, front.cutoff = 0) {
  
  data.df <- word.freq(s)
  
  # End Cutoff
  data.df <- data.df[data.df$Rank <= end.cutoff,]
  
  # Front Cutoff
  data.df <- data.df[data.df$Rank > front.cutoff,]
  data.df$Rank <- data.df$Rank - front.cutoff
  
  # Normalize
  data.df$Freq.n <- (data.df$Freq - min(data.df$Freq)) / (max(data.df$Freq) - min(data.df$Freq))

  # Rename columns to pick 'x' and 'y'
  colnames(data.df) <- c('Word', 'Freq', 'Prop', 'x', 'y')

  # Starting parameters
  start <- list(beta = -1)

  # Model
  model <- nls(y ~ x ^ beta, data = data.df, start = start)

  model.sum <- summary(model)
  beta <- model.sum$coefficients[1]
  beta.error <- model.sum$coefficients[2]
  
  df <- data.frame(beta, beta.error)

  return (df) 
  
}









# This is a function which runs all of the statistics that I think could
# be relevant. It is used to create the statistics files. Update the 
# stats here.


multi_stats <- function (doc) {
  
  # Word Properties
  doc.list <- unlist(strsplit(doc, " "))
  word.len <- length(doc.list)
  doc.table <- table(doc.list)
  doc.table <- doc.table / sum(doc.table)
  doc.table <- doc.table[order(doc.table, decreasing = TRUE)]
  word.size <- length(doc.table)
  
  # Character Properties
  char.len <- str_length(doc)
  char.size <- length(table(unlist(strsplit(doc, ""))))
  
  # Character Entropy
  char.entropy <- sumentropy(doc)
  char.h0 <- as.numeric(char.entropy[2])
  char.h1 <- as.numeric(char.entropy[3])
  char.H2 <- as.numeric(char.entropy[4])
  char.h2 <- as.numeric(char.entropy[5])
  char.H3 <- as.numeric(char.entropy[8])
  char.h3 <- as.numeric(char.entropy[9])
  
  # Word Entropy
    # MAWE 10000
  word.mawe10000 <- mawe(doc, window=10000)
  word.mawe10000.h0 <- as.numeric(word.mawe10000[2])
  word.mawe10000.h1 <- as.numeric(word.mawe10000[3])
  word.mawe10000.H2 <- as.numeric(word.mawe10000[4])
  word.mawe10000.h2 <- as.numeric(word.mawe10000[5])
  word.mawe10000.H3 <- as.numeric(word.mawe10000[6])
  word.mawe10000.h3 <- as.numeric(word.mawe10000[7])
  
  
  
  # MATTR
  mattr500 <- as.numeric(mattr(doc, window = 500))
  mattr1000 <- as.numeric(mattr(doc, window = 1000))
  mattr2000 <- as.numeric(mattr(doc, window = 2000))
  
  # Frequent Words
  word1 <- names(doc.table[1])
  word2 <- names(doc.table[2])
  word3 <- names(doc.table[3])
  word4 <- names(doc.table[4])
  word5 <- names(doc.table[5])
  
  # Word Proportions
  word1.prop <- sum(doc.table[1])
  word2.prop <- sum(doc.table[2])
  word5.prop <- sum(doc.table[1:5])
  word10.prop <- sum(doc.table[1:10])
  word25.prop <- sum(doc.table[1:25])
  word50.prop <- sum(doc.table[1:50])
  word100.prop <- sum(doc.table[1:100])
  
  # Zipf Fitting
  doc.zipf <- zipf_fit(doc, end.cutoff = 1000000000)
  beta <- doc.zipf$beta
  beta.error <- doc.zipf$beta.error
  
  doc.zipf.100 <- zipf_fit(doc, end.cutoff = 100)
  beta.100 <- doc.zipf.100$beta
  beta.100.error <- doc.zipf.100$beta.error
  
  
  df <- data.frame(word.len, word.size, char.len, char.size, char.h0, char.h1, char.H2, char.h2, char.H3, char.h3, word.mawe10000.h0, word.mawe10000.h1, word.mawe10000.H2, word.mawe10000.h2, word.mawe10000.H3, word.mawe10000.h3, mattr500, mattr1000, mattr2000, word1, word2, word3, word4, word5, word1.prop, word2.prop, word5.prop, word10.prop, word25.prop, word50.prop, word100.prop, beta, beta.error, beta.100, beta.100.error)
  
  return (df)  
}




# MAWE Moving Average Word Entropy: 
# average word entropy statistics within a moving window

mawe <- function(doc, window = 500) {
  
  word.list <- unlist(strsplit(doc, " "))
  
  # Super-useful partition function for R:
  word.part <- split(word.list, as.integer((seq_along(word.list) - 1) / window))
  
  # Type-token frequency for each partion
  word.num.list <- c()
  h0.list <- c()
  h1.list <- c()
  H2.list <- c()
  h2.list <- c()
  H3.list <- c()
  h3.list <- c()
  
  for (x in word.part) {
    
    x.whole <- concatenate(x)
    x.freq = table(x)
    word.num = length(x.freq)
    x.norm <- x.freq/sum(x.freq)
    
    x.freq.bi.ngram = ngram(x.whole, n=2)
    x.freq.bi = get.phrasetable(x.freq.bi.ngram)
    x.norm.bi <- x.freq.bi$prop
    
    x.freq.tri.ngram = ngram(x.whole, n=3)
    x.freq.tri = get.phrasetable(x.freq.tri.ngram)
    x.norm.tri <- x.freq.tri$prop
    
    h0 = log2(length(x.freq))
    
    h1 = -sum(log2(x.norm)*x.norm)
    
    H2 = -sum(log2(x.norm.bi)*x.norm.bi)
    
    H3 = -sum(log2(x.norm.tri)*x.norm.tri)
    
    h2 = H2-h1
    
    h3 = H3-H2
    
    word.num.list <- c(word.num.list, word.num)
    h0.list <- c(h0.list, h0)
    h1.list <- c(h1.list, h1)
    H2.list <- c(H2.list, H2)
    h2.list <- c(h2.list, h2)
    H3.list <- c(H3.list, H3)
    h3.list <- c(h3.list, h3)    
    
  }
  
  # Return the avg
  s.summary = as.table(c(mean(word.num.list), mean(h0.list), mean(h1.list), mean(H2.list), mean(h2.list), mean(H3.list), mean(h3.list)))
  rownames(s.summary) = c("word.num","h0", "h1", "H2", "h2", "H3", "h3")
  s.summary
}



