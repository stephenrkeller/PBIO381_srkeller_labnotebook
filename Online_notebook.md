# Project title    

### Author: Stephen Keller     
### Ecological Genomics 2017:   

## Overall Description of notebook (fill in with description)  

## Date started: (2017-01-17)   
## Date end:   (year-month-day)    

## Philosophy   
Science should be reproducible and one of the best ways to achieve this is by logging research activities in a notebook. Because science/biology has increasingly become computational, it is easier to document computational projects in an electronic form, which can be shared online through Github.    

### Helpful features of the notebook     

It is absolutely critical for your future self and others to follow your work. The notebook is set up with a series of internal links from the table of contents. All notebooks should have a table of contents which has the "Page", date, and title (information that allows the reader to understand your work). Also, one of the perks of keeping all activities in a single document is that you can search and find elements quickly. You can document anything you'd like, aside from logging your research activities. For example, feel free to log all/any ideas for your research project([example](https://github.com/adnguyen/Notebooks_and_Protocols/blob/master/2016_notebook.md#page-39-2016-06-13-post-doc-project-idea-assessing-current-impacts-of-climate-change-in-natural-populations)) as an entry, or write down notes for a paper([example](https://github.com/adnguyen/Notebooks_and_Protocols/blob/master/2016_notebook.md#id-section36). Lastly, you can share specific entries because of the three "#" automatically creates a link when the notebook renders on github.      


### Table of contents for 60 entries (Format is *Page: Date(with year-month-day). Title*)        
* [Page 1: 2017-01-17](#id-section1). Initial entry
* [Page 2:](#id-section2).
* [Page 3:](#id-section3).
* [Page 4:](#id-section4).
* [Page 5:](#id-section5).
* [Page 6:](#id-section6).
* [Page 7:](#id-section7).
* [Page 8:](#id-section8).
* [Page 9:](#id-section9).
* [Page 10:](#id-section10).
* [Page 11:](#id-section11).
* [Page 12:](#id-section12).
* [Page 13:](#id-section13).
* [Page 14:](#id-section14).
* [Page 15:](#id-section15).
* [Page 16:](#id-section16).
* [Page 17:](#id-section17).
* [Page 18:](#id-section18).
* [Page 19:](#id-section19).
* [Page 20:](#id-section20).
* [Page 21:](#id-section21).
* [Page 22:](#id-section22).
* [Page 23:](#id-section23).
* [Page 24:](#id-section24).
* [Page 25:](#id-section25).
* [Page 26:](#id-section26).
* [Page 27:](#id-section27).
* [Page 28:](#id-section28).
* [Page 29:](#id-section29).
* [Page 30:](#id-section30).
* [Page 31:](#id-section31).
* [Page 32:](#id-section32).
* [Page 33:](#id-section33).
* [Page 34:](#id-section34).
* [Page 35:](#id-section35).
* [Page 36:](#id-section36).
* [Page 37:](#id-section37).
* [Page 38:](#id-section38).
* [Page 39:](#id-section39).
* [Page 40:](#id-section40).
* [Page 41:](#id-section41).
* [Page 42:](#id-section42).
* [Page 43:](#id-section43).
* [Page 44:](#id-section44).
* [Page 45:](#id-section45).
* [Page 46:](#id-section46).
* [Page 47:](#id-section47).
* [Page 48:](#id-section48).
* [Page 49:](#id-section49).
* [Page 50:](#id-section50).
* [Page 51:](#id-section51).
* [Page 52:](#id-section52).
* [Page 53:](#id-section53).
* [Page 54:](#id-section54).
* [Page 55:](#id-section55).
* [Page 56:](#id-section56).
* [Page 57:](#id-section57).
* [Page 58:](#id-section58).
* [Page 59:](#id-section59).
* [Page 60:](#id-section60).

------
<div id='id-section1'/>
### Page 1: 2017-01-17.  Initial entry   


To embed an image, start with the "!", then [], then the URL in ()
```
![](URL)
```
![](https://cloud.githubusercontent.com/assets/12184909/22031778/c2232e86-dcaf-11e6-90c9-782206dcd824.png)

This image is a fastqc plot of Brittany's spruce GBS data

Here's an example of pasting output from the commandline:

```
stephen@Wright:.../Centaurea_GBS/fastq$ zcat srkeller_GBS_20160812_C1_keller_R1_001.fastq.gz | head
@HISEQ:140:160817_SNL128_0140_AC96CYACXX:8:1101:1425:1938 1:N:0:
NTTGTGACAGCTTCCTGAAAATCAAACAAAGAAATTAACAAACCGATCATGAGATCGGTAGACTTAATTCATGAGAGATTTACCTGAAAAAACGTGGCAGA
+
#0<<BBFFFFFFFFFFFFBFBFFFBBBBF0F<FBFIFFFFFBFFF7BFFIIIBBBFBFBFFFF<BFFFBFBFFBBBBBBBBBFBBBBBB<B<BBBBBBBBB
@HISEQ:140:160817_SNL128_0140_AC96CYACXX:8:1101:1398:1945 1:N:0:
NATGTCAGCGCTTCATGAAGCGCACCGGCAAGATGCCCAACATGATCCACGTCGGCACCTATTCGGCGGTGCTGAGATCGGAAGAGCGGTTCAGCAGGAAT
+
#0<FFFFFFFFFFIIIIIIIIIIIIIIIIIFFIFIIIIIFFIFIIIIFFFFFBFFFBFFBFFFFFFFFF0<BFF<7BBBFFFBBFBBFF77BBFBBFFFF<
@HISEQ:140:160817_SNL128_0140_AC96CYACXX:8:1101:1292:1959 1:N:0:
NCGGTAATACAGCGCCGACATGGAGCCCCTGCTGTAGCCATGGCTTGCTGGTTGCCCTCAATCCGGCCTAGTCAGCTGAGATCGGAAGAGCGGTTCAGCAG
stephen@Wright:.../Centaurea_GBS/fastq$ ll
total 12991974
drwxrwx--- 2 stephen klab           0 Jan  6 15:32 ./
drwxrwx--- 2 stephen klab           0 Jan 12 15:00 ../
drwxrwx--- 2 stephen klab           0 Jan  6 13:43 distribs/
-rw-rw---- 1 stephen klab        4096 Dec 22 12:36 ._.DS_Store
-rw-rw---- 1 stephen klab       10244 Jan 12 16:57 .DS_Store
drwxrwx--- 2 stephen klab           0 Jan  6 15:32 fastqc/
drwxrwx--- 2 stephen klab           0 Jan  5 17:33 old/
drwxrwx--- 2 stephen klab           0 Jan 12 16:57 parsed/
-rw-rw---- 1 stephen klab 13303765656 Jan  5 21:52 srkeller_GBS_20160812_C1_keller_R1_001.fastq.gz
drwxrwx--- 2 stephen klab           0 Jan  6 13:41 summaries/
stephen@Wright:.../Centaurea_GBS/fastq$ 

```


------
<div id='id-section2'/>
### Page 2: 
------
<div id='id-section3'/>
### Page 3:
------
<div id='id-section4'/>
### Page 4:
------
<div id='id-section5'/>
### Page 5:
------
<div id='id-section6'/>
### Page 6:
------
<div id='id-section7'/>
### Page 7:
------
<div id='id-section8'/>
### Page 8:
------
<div id='id-section9'/>
### Page 9:
------
<div id='id-section10'/>
### Page 10:
------
<div id='id-section11'/>
### Page 11:
------
<div id='id-section12'/>
### Page 12:
------
<div id='id-section13'/>
### Page 13:
------
<div id='id-section14'/>
### Page 14:
------
<div id='id-section15'/>
### Page 15:
------
<div id='id-section16'/>
### Page 16:
------
<div id='id-section17'/>
### Page 17:
------
<div id='id-section18'/>
### Page 18:
------
<div id='id-section19'/>
### Page 19:
------
<div id='id-section20'/>
### Page 20:
------
<div id='id-section21'/>
### Page 21:
------
<div id='id-section22'/>
### Page 22:
------
<div id='id-section23'/>
### Page 23:
------
<div id='id-section24'/>
### Page 24:
------
<div id='id-section25'/>
### Page 25:
------
<div id='id-section26'/>
### Page 26:
------
<div id='id-section27'/>
### Page 27:
------
<div id='id-section28'/>
### Page 28:
------
<div id='id-section29'/>
### Page 29:
------
<div id='id-section30'/>
### Page 30:
------
<div id='id-section31'/>
### Page 31:
------
<div id='id-section32'/>
### Page 32:
------
<div id='id-section33'/>
### Page 33:
------
<div id='id-section34'/>
### Page 34:
------
<div id='id-section35'/>
### Page 35:
------
<div id='id-section36'/>
### Page 36:
------
<div id='id-section37'/>
### Page 37:
------
<div id='id-section38'/>
### Page 38:
------
<div id='id-section39'/>
### Page 39:
------
<div id='id-section40'/>
### Page 40:
------
<div id='id-section41'/>
### Page 41:
------
<div id='id-section42'/>
### Page 42:
------
<div id='id-section43'/>
### Page 43:
------
<div id='id-section44'/>
### Page 44:
------
<div id='id-section45'/>
### Page 45:
------
<div id='id-section46'/>
### Page 46:
------
<div id='id-section47'/>
### Page 47:
------
<div id='id-section48'/>
### Page 48:
------
<div id='id-section49'/>
### Page 49:
------
<div id='id-section50'/>
### Page 50:
------
<div id='id-section51'/>
### Page 51:
------
<div id='id-section52'/>
### Page 52:
------
<div id='id-section53'/>
### Page 53:
------
<div id='id-section54'/>
### Page 54:
------
<div id='id-section55'/>
### Page 55:
------
<div id='id-section56'/>
### Page 56:
------
<div id='id-section57'/>
### Page 57:
------
<div id='id-section58'/>
### Page 58:
------
<div id='id-section59'/>
### Page 59:
------
<div id='id-section60'/>
### Page 60:

------
