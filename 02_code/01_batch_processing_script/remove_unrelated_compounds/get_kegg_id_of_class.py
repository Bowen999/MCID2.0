import requests
from bs4 import BeautifulSoup
import re

# url = "https://www.genome.jp/brite/br08002"
# response = requests.get(url)
# response.raise_for_status()  # Ensure the request was successful

# # Create a BeautifulSoup object
# soup = BeautifulSoup(response.text, 'html.parser')


import re

# Load the file content
with open('/Users/bowen/Desktop/soup.txt', 'r', encoding='utf-8') as file:
    content = file.read()

# Regular expression to find the IDs
pattern = re.compile(r'\bC\d{5}\b')

# Find all occurrences of the pattern
entry_ids = re.findall(pattern, content)

# Print the list of entry IDs
print(entry_ids)


