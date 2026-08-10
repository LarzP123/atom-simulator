cabal test --enable-coverage

$htmlDir = Get-ChildItem -Path "dist-newstyle" -Recurse -Filter "html" -Directory | Select-Object -First 1
if ($htmlDir) {
    Remove-Item testdata -Recurse -Force -ErrorAction SilentlyContinue
    New-Item testdata -ItemType Directory -Force | Out-Null
    Get-ChildItem $htmlDir.FullName |
        Where-Object { $_.Name -notlike "hpc_index_*" } |
        Copy-Item -Destination testdata -Recurse
    
    # Replace hpc_index_*.html references in hpc_index.html
    $hpcIndexPath = Join-Path testdata "hpc_index.html"
    if (Test-Path $hpcIndexPath) {
        $content = Get-Content $hpcIndexPath -Raw
        $content = $content -replace 'hpc_index_[^"''>\s]+\.html', 'hpc_index.html'
        Set-Content $hpcIndexPath -Value $content
        Write-Host "Updated hpc_index.html references"
    }
    Write-Host "Coverage report at $(pwd)\testdata"
} else {
    Write-Host "No coverage found"
}