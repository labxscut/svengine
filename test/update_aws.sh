#configure AWS_PATH if needed
#export AWS_PATH="/home/bkd/.local/bin/"
#export PATH=$AWS_PATH:$PATH

#update svengine ami as needed
ssh bkd@ec2-54-202-133-252.us-west-2.compute.amazonaws.com

#now clone back to sake and upload a new svengine test package
tar -czvf svengine.test.full.tgz -T svengine.test.full.list
aws s3 cp svengine.test.full.tgz s3://lixiabucket/svengine.test.full.tgz --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers 
#need to change permission every time it is updated

#now redeposit a new ova
#aws ec2 create-instance-export-task --instance-id i-0bcc981a6b1ab6744 --target-environment vmware --export-to-s3-task DiskImageFormat=vmdk,ContainerFormat=ova,S3Bucket=lixiabucket,S3Prefix=svengineDemo.2017Apr26.
taskid=$(aws ec2 create-instance-export-task --instance-id i-04923cf01efe71033 --target-environment vmware --export-to-s3-task DiskImageFormat=vmdk,ContainerFormat=ova,S3Bucket=lixiabucket,S3Prefix=svengineDemo. 2>&1 | grep -m1 -o export-i-[0-9a-z]*)
echo $taskid
aws ec2 describe-export-tasks --export-task-ids $taskid #show status
aws s3 cp s3://lixiabucket/svengineDemo.$taskid.ova s3://lixiabucket/svengineDemo.ova --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
#need to change permission every time it is updated
