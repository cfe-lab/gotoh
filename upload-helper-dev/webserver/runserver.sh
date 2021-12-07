#!/bin/bash

# run the web server in development mode
uvicorn webservermain:app --reload --host 0.0.0.0
