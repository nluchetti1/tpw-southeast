name: Update Weather Dashboard
on:
  schedule:
    - cron: '*/15 * * * *'
  workflow_dispatch:

permissions:
  contents: write

jobs:
  run-analysis:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout Repo
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.10'

      - name: Install System Geo-Libs
        run: |
          sudo apt-get update
          # Install essential engines for GRIB2 (RAP) and Mapping (Cartopy)
          sudo apt-get install -y libgeos-dev libproj-dev proj-data proj-bin libhdf5-dev libeccodes-dev libeccodes-tools

      - name: Install Python Dependencies
        run: |
          python -m pip install --upgrade pip
          # Force a clean install of Herbie and its GRIB2 backend (cfgrib)
          pip install --no-cache-dir herbie cfgrib satpy pyresample xarray netCDF4 metpy cartopy siphon matplotlib scipy requests

      - name: Diagnostic Check
        run: |
          echo "Checking for Herbie..."
          python -c "import herbie; print('SUCCESS: Herbie is installed at', herbie.__file__)"
          echo "Checking for GRIB engine..."
          python -c "import cfgrib; print('SUCCESS: cfgrib is ready')"

      - name: Generate Plot
        run: python generate_map.py

      - name: Commit Update
        run: |
          git config --local user.email "action@github.com"
          git config --local user.name "GitHub Action"
          git add output_map.png
          git commit -m "Auto-Update: $(date)" || exit 0
          git push
