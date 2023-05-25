# Example adapter sequence file

# Example metadata file

The metadata file is a nested dictionary in json format. The outer dictionary contains library names and genome as dictionary keys. The value for each library name is a dictionary whose **key**-**value** pairs are:

1. **key**: prefix; **value**: the filename prefix assigned to this library as given in the filenames,usually the library name itself. 
2. **key**: reps; **value**: all biological and technical replicates of the library as provided in the filenames. The bio and tech replicate names should be separated by an underscore. Unique bio and tech rep combination should be separated by a space.  
3. **key**: readpairs; **value**: the read pair names as given in the filename usually R1 and R2. Paired end read names should be separated by space.
4. **key**: fastqsuffix; **value**: the filename suffix usually "fastq.gz" or "fq.gz".
5. **key**: shortform; **value**: the library shortform to be assigned by the user based on which new subdirectories to store results will be created.

The value for the genome key is a dictionary whose **key**-**value** pairs are:

1. **key**: fasta; **value**: the filepath of the genome fasta file.
2. **key**: gtf; **value**: the filepath of the gtf file.
