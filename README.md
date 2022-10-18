# t3dbms2:Tool for manipulating T3DB databases
# Installation
You can install t3dbms2 from Github

if(!require(remotes)){
install.packages("remotes")
}
remotes::install_github("Zhujiamin/t3dbms2")

# Usage
get_t3db(page_num,name,sleep_time)

page_num:The number of pages you want to crawl
name:The class of substance, such as "pollutant"
sleep_time:The time between each acquisition of data

# other
There will be problems, improvement



