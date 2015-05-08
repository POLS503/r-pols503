#' Convert column to dummy variables
#'
#' \code{to_dummies} converts a column into multiple columns of dummy variables.
#' \code{to_dummies_} is the standard evaluation version of this function.
#'
#' @param .data A data frame
#' @param col The bare column name
#' @param remove If \code{TRUE}, then remove the input column from the output
#'   data frame.
#' @param drop_level If \code{TRUE}, then drop a level when creating the variables.
#'   This is the default for models.
#' @param vars A \code{character} vector with the names of the new variables,
#'   or a function with two arguments, the column name and the factor levels, which
#'   returns a character vector of new variable names.
#' @return A data frame.
#' @examples
#' example_df <- data.frame(x = c(rep("a", 2), rep("b", 2)))
#' to_dummies(example_df, x)
#' to_dummies(example_df, x, remove = FALSE)
#' to_dummies(example_df, x, drop_level = TRUE)
#' @export
to_dummies <- function(.data, col,
                       remove = TRUE, drop_level = FALSE,
                       vars = NULL) {
  col <- col_name(substitute(col))
  to_dummies_(.data, col, remove = remove,
              drop_level = drop_level)
}

#' @rdname to_dummies
#' @export
to_dummies_ <- function(.data, col,
                        remove = TRUE, drop_level = FALSE,
                        vars = NULL) {
  UseMethod("to_dummies_")
}

to_dummies_.tbl_df <- function(.data, col,
                               remove = TRUE, drop_level = FALSE,
                               vars = NULL) {
  dplyr::tbl_df(NextMethod())
}

has_attr <- function(x, which) {
  ! is.null(attr(x, which))
}

fill_na <- function(x, fill = 0) {
  x[is.na(x)] <- fill
  x
}

to_dummies_.data.frame <- function(.data, col,
                                   drop_level = FALSE,
                                   remove = TRUE,
                                   vars = NULL,
                                   contr = NULL,
                                   na_level = FALSE,
                                   na_varname = NULL) {
  stopifnot(is.character(col), length(col) == 1)
  # do this so that the column keeps any factor attributes
  if (! is.factor(.data[[col]])) {
    x <- as.factor(as.character(.data[[col]]))
  } else {
    x <- .data[[col]]
  }
  xlvls <- levels(x)
  if (is.null(contr)) {
    if (! has_attr(x, "contrasts")) {
      contr <- contr.treatment(xlvls, contrasts = drop_level)
    } else {
      contr <- attr(x, "contrasts")
    }
  }
  if (is.null(rownames(contr))) {
    rownames(contr) <- xlvls
  }
  if (!is.null(vars)) {
    if (is.function(vars)) {
      var_names <- vars(col, xlvls)
    } else {
      var_names <- as.character(vars)
    }
  } else {
    var_names <- paste0(col, xlvls)
  }
  if (! na_level) {
    if (!is.null(na_varname)) {
      if (is.function(na_varname)) {
        na_x <- na_varname(col)
      } else {
        na_x <- as.character(na_varname)
      }
    } else {
      na_x <- paste0(col, NA)
    }
  }
  for (i in seq_len(ncol(contr))) {
    newname <- var_names[i]
    .data[[newname]] <- contr[as.character(x), i]
  }
  if (remove) {
    .data[[col]] <- NULL
  }
  .data
}

which_w_default <- function(x, y, default = NA) {
  ret <- x[! is.na(y) & as.logical(y)]
  if (length(ret) == 0) ret <- default
  else ret <- ret[1]
  ret
}

#' Convert dummy variable columns to a categorical variable column
#'
#' \code{from_dummies} converts dummy variable columns to a categorical variable
#' column.
#' \code{from_dummies_} is the standard evaluation version.
#'
#' @param .data A data frame
#' @param col The bare column name which will be created.
#' @param ...  Comma separated list of unquoted expressions.
#'   You can treat variable names like they are positions.
#'   Use positive values to select variables; use negative values to drop variables.
#' @param .dots	Use \code{from_dummies_()} to do standard evaluation.
#' @param remove If \code{TRUE}, then remove the input column from the output
#'   data frame.
#' @param default The value to use if no dummy variables match.
#' @return A data frame.
#' @examples
#' example_df <- data.frame(a = c(1, 0, 0), b = c(0, 1, 0))
#' from_dummies(example_df, x, a, b)
#' from_dummies(example_df, x, a, b, remove = FALSE)
#' from_dummies(example_df, x, a, b, default = "c")
#' @rdname from_dummies
#' @export
from_dummies <- function(.data, col, ..., default = NA_character_,
                         remove = TRUE) {
  col <- col_name(substitute(col))
  from_dummies_(.data, col, .dots = lazyeval::lazy_dots(...),
                default = default, remove = remove)
}

#' @rdname from_dummies
#' @export
from_dummies_ <- function(.data, col, ..., .dots, default = NA_character_,
                          remove = TRUE) {
  UseMethod("from_dummies_")
}

from_dummies_.data.frame <- function(.data, col, ..., .dots,
                                     default = NA_character_, remove = TRUE) {
  stopifnot(is.character(col), length(col) == 1)
  dots <- lazyeval::all_dots(.dots, ...)
  from <- select_vars_(names(.data), dots)
  catvar <-
    apply(.data[ , from, drop = FALSE], 1,
          function(x) which_w_default(from, x, default = default))
  .data <- append_col(.data, catvar, col)
  if (remove) {
    for (i in from) {
      .data[[i]] <- NULL
    }
  }
  .data
}

from_dummies_.tbl_df <- function(.data, col, from, default = NA_character_,
                                 remove = TRUE) {
  dplyr::tbl_df(NextMethod())
}


#' Transform a data frame using a formula
#'
#' This uses \code{model.matrix} formula to transform a data frame.
#'
#' @param .data A data frame
#' @param formula A formula
#' @param append Add columns to \code{.data}. If \code{FALSE}, then return a new
#'   data frame.
#' @param na.action The \code{\link{na.action}} to use when creating the \code{model.matrix}.
#' @return A data frame with the transformations specified in the formula.
#' @export
model_df <- function(.data, formula, append = FALSE, na.action = "na.pass") {
  orig_na_action <- getOption("na.action")
  on.exit(options(na.action = orig_na_action))
  options(na.action = na.action)
  df <- as.data.frame(model.matrix(formula, data = .data))
  if (append) {
    for (i in names(df)) {
      .data[[i]] <- df[[i]]
    }
    df <- .data
  }
  df
}
