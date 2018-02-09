1) Accession ID is presented management command
2) Management command creates SurveyorJob
3) SurveyorJob reads API endpoint for that accession ID
4) SurveyorJob populates Experiment
5) Surveyor Job creates DownloaderJobs. Each downloader is told the number of expected samples.
6) Each DownloaderJob downloads a sample file. 
    If an individual file, creates a Sample item and appropriate Sample Association.
    If a zip, unpacks and creates multiple Sample items and creates appropriate ProcessorJobs.
7) ProcessorJobs execute. When they finish, ComputationalResults are created.
    If the number of outputs in the directory is correct / there is a lock file, create multi-Sample processor jobs.
    (If this is a final processor job, uploader jobs are created. TODO this.)
    