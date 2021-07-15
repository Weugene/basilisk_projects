#!/bin/bash

# 0. Ensure AWS CLI tool installed: pip install awscli
# 1. Make shell script executable: chmod u+x ec2-spot-prices.sh
# 2. Run script and provide instance type to check: ./ec2-spot-prices.sh c4.8xlarge
# 3. Script runs and outputs full stops while querying the Amazon API,
#      returning three columns: Region+AZ, Instance Type, Current Spot Price in $

allSpot=""
for Reg in eu-west-1 eu-west-2 eu-west-3 eu-central-1 eu-north-1 ca-central-1 us-east-1 us-east-2 us-west-1 us-west-2 sa-east-1 ap-southeast-1 ap-northeast-1 ap-east-1 ap-northeast-2 ap-southeast-2 ap-south-1
do
  printf "."
  nextSpot=$(aws --region $Reg ec2 describe-spot-price-history \
               --filters '[{"Name":"product-description","Values":["Linux/UNIX (Amazon VPC)"]}]' \
               --instance-types $@ \
               --start-time `date "+%Y-%m-%dT%H:%M:00"` \
               --output text \
             | cut -f 2,3,5)
  allSpot="$allSpot\n$nextSpot"
done;
printf "\n"
printf "$allSpot" | sort -gk 3,3 | grep -v '^$'
