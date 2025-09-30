#make sure the ecr repository exists and awscli is installed and configured with the keys


aws ecr get-login-password --region ap-southeast-1 | docker login --username AWS --password-stdin 211125343458.dkr.ecr.ap-southeast-1.amazonaw
s.com

docker tag mi-heart-st:latest 211125343458.dkr.ecr.ap-southeast-1.amazonaws.com/mi-heart-st:latest

docker push 211125343458.dkr.ecr.ap-southeast-1.amazonaws.com/mi-heart-st:latest