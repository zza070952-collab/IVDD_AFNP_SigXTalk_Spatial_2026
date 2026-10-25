param([Parameter(Mandatory=$true)][string]$ProjectRoot,[Parameter(Mandatory=$true)][string]$Docx,[Parameter(Mandatory=$true)][string]$Pdf,[Parameter(Mandatory=$true)][string]$Receipt)
$ErrorActionPreference = 'Stop'
$inputDoc = (Resolve-Path -LiteralPath $Docx).Path
$targetPdf = [IO.Path]::GetFullPath($Pdf)
$recordPath = [IO.Path]::GetFullPath($Receipt)
if ($inputDoc -eq $targetPdf -or $inputDoc -eq $recordPath -or $targetPdf -eq $recordPath) { throw 'Input, PDF and receipt must be separate files' }
if ((Test-Path -LiteralPath $targetPdf) -or (Test-Path -LiteralPath $recordPath)) { throw 'Choose fresh PDF and receipt destinations' }
$taskRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path + [IO.Path]::DirectorySeparatorChar
$checkedDirectory = [IO.DirectoryInfo]::new((Split-Path -Parent $PSCommandPath))
while ($null -ne $checkedDirectory) {
 if (Test-Path -LiteralPath (Join-Path $checkedDirectory.FullName '.git')) {
  $codeRoot = $checkedDirectory.FullName + [IO.Path]::DirectorySeparatorChar
  if ($taskRoot.StartsWith($codeRoot,[StringComparison]::OrdinalIgnoreCase)) { throw 'Private document outputs must be outside this publication checkout' }
  break
 }
 $checkedDirectory = $checkedDirectory.Parent
}
foreach ($candidate in @($inputDoc,$targetPdf,$recordPath)) {
 if (-not $candidate.StartsWith($taskRoot,[StringComparison]::OrdinalIgnoreCase)) { throw 'Every document/export target must stay in the explicit private project root' }
}
$before = (Get-FileHash -LiteralPath $inputDoc -Algorithm SHA256).Hash.ToLower()
[IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($targetPdf)) | Out-Null
[IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($recordPath)) | Out-Null
$wordApp = New-Object -ComObject Word.Application
$wordApp.Visible = $false
$wordApp.DisplayAlerts = 0
$document = $null
try {
 $document = $wordApp.Documents.Open($inputDoc,$false,$true,$false)
 $document.Repaginate()
 $pages = $document.ComputeStatistics(2)
 $document.ExportAsFixedFormat($targetPdf,17)
 $wordVersion = $wordApp.Version
} finally {
 if ($null -ne $document) { $document.Close(0); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($document) }
 $wordApp.Quit(0)
 [void][Runtime.InteropServices.Marshal]::ReleaseComObject($wordApp)
}
$after = (Get-FileHash -LiteralPath $inputDoc -Algorithm SHA256).Hash.ToLower()
if ($before -ne $after) { throw 'Input DOCX changed during read-only PDF export' }
if (-not (Test-Path -LiteralPath $targetPdf)) { throw 'Word did not export a PDF' }
$result = [ordered]@{status='EXPORTED_PENDING_PAGE_VISUAL_QA';exporter='Microsoft Word native read-only ExportAsFixedFormat';word_version=$wordVersion;pages=$pages;docx_sha256=$after;pdf_sha256=(Get-FileHash -LiteralPath $targetPdf -Algorithm SHA256).Hash.ToLower();docx_unchanged=$true;user_Word_sessions_closed=$false}
$result | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $recordPath -Encoding utf8
$result | ConvertTo-Json -Compress
