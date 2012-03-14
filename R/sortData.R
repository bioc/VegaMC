sortData <- function(dataset, output_file_name = "")
{
    if (output_file_name == "") {
        message("ERROR: You must provide the name of the file",
                "in which sorted data are saved")
    } else {
        res <- .C("sort_data", as.character(dataset),
                  as.character(output_file_name))
    }
}
