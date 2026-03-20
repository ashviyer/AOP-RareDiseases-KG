import requests

API_KEY = "YOUR_API_KEY"
HEADERS = {
    "Authorization": API_KEY,
    "accept": "application/json"
}

UMLS_CUI = "C0020179"  # Huntington's disease

url = "https://api.disgenet.com/api/v1/gda/summary"

response = requests.get(
    url,
    params={"disease_id": UMLS_CUI},
    headers=HEADERS,
    timeout=15
)

data = response.json()

print("STATUS:", data.get("status"))

# Handle success vs error safely
if data.get("status") == "OK":
    paging = data.get("paging", {})
    payload = data.get("payload", [])

    print("TOTAL RESULTS:", paging.get("totalElements", 0))
    print("RESULTS IN PAGE:", paging.get("totalElementsInPage", 0))

    print("\nFirst result example:")
    if payload:
        print(payload[0])
else:
    print("ERROR DETAILS:", data.get("payload", {}).get("details"))