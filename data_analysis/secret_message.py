import requests
from bs4 import BeautifulSoup

def get_message_from_url(url):
    # Fetch the document content
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")

    # Locate the table containing the data
    tables = soup.find_all("table")

    # Extract data from the table (assumes the first table is the relevant one)
    data = []
    for table in tables:
        rows = table.find_all("tr")
        for row in rows[1:]:  # Skip the header row
            columns = row.find_all(["td", "th"])
            if len(columns) >= 3:  # Ensure we have enough columns
                x = int(columns[0].get_text(strip=True))
                char = columns[1].get_text(strip=True)
                y = int(columns[2].get_text(strip=True))
                data.append({"x": x, "y": y, "char": char})

    # Determine grid size
    max_x = max(item["x"] for item in data) + 1
    max_y = max(item["y"] for item in data) + 1

    # Create an empty grid
    grid = [[" " for _ in range(max_x)] for _ in range(max_y)]

    # Populate the grid with reversed y-coordinate
    for item in data:
        x, y, char = item["x"], item["y"], item["char"]
        grid[max_y - 1 - y][x] = char  # Flip the y-coordinate

    # Print the grid
    # Convert the grid to a single string variable
    grid_output = "\n".join("".join(row) for row in grid)

    return print(grid_output)

# Print the grid (optional, for verification)


url = "https://docs.google.com/document/d/e/2PACX-1vQGUck9HIFCyezsrBSnmENk5ieJuYwpt7YHYEzeNJkIb9OSDdx-ov2nRNReKQyey-cwJOoEKUhLmN9z/pub"

get_message_from_url(url)
