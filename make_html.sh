#!/bin/bash

parallel -uj 4 md_to_html {} ::: *.md
