setClass("Snakemake", representation(
    input = "list",
    output = "list",
    params = "list",
    wildcards = "list",
    log = "list",
    config = "list"
))

make_list <- function(item) {
    if (typeof(item) == "list") {
        return(item)
    }
    return(list(item))
}

snakemake_debug <- function(
                            input = list(),
                            output = list(),
                            params = list(),
                            wildcards = list(),
                            log = list(),
                            config = list()) {
    return(new("Snakemake",
        input = make_list(input),
        output = make_list(output),
        params = make_list(params),
        wildcards = make_list(wildcards),
        log = make_list(log),
        config = make_list(config)
    ))
}
