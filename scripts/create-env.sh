#!/bin/bash

echo "# Created by running \"bash ./scripts/create-env.sh\"
USER_ID=$(id -u $USER)
USER=$USER
" >| .env

echo ".env written to .env"