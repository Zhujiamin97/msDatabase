#                   生活不止眼前的苟且,还有吃和远方

![image](https://user-images.githubusercontent.com/93595586/196358703-6e844094-03e0-44db-803a-5e039d6b761c.png)


# msdatabase:Tool for manipulating Open source mass spectrometry databases
# Installation
You can install msdatabase from Github

if(!require(remotes)){
install.packages("remotes")
}

remotes::install_github("Zhujiamin/msdatabase")

# Usage
get_t3db(page_num,name,sleep_time)

start_page_num:The page number you want to start crawling

end_page_num:The page number you want to end the crawl

name:The class of substance, such as "pollutant"

Please select a category from the list below

![image](https://user-images.githubusercontent.com/93595586/196375621-2c955dad-aa7a-463e-b4d9-9b87fed3e442.png)

sleep_time:The time between each acquisition of data

![image](https://user-images.githubusercontent.com/93595586/196357074-4fee4e08-b667-451f-9cb5-5b90cbc4cc04.png)

# other tools HMDB
Use the get_hmdb_id(file) function to obtain the HMDBID corresponding to the metabolite in the xml file

Use hmdb_ms2 function to obtain the corresponding spectrum of metabolites

# warn
It is best to crawl multiple times, otherwise access will be rejected by the website

There will be problems, improvement



