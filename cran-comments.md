## Resubmission
This is a resubmission.

## NOTES
> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Done. We added the main reference of the package in the description field of the
DESCRIPTION file.

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      envelope.Rd: \value
      methodsPLreg.Rd: \value
      plot.PLreg.Rd: \value
      residuals.PLreg.Rd: \value

Done. We added \\value to envelope.Rd, plot.PLreg.Rd and residuals.PLreg.Rd. 
For the methodsPLreg.Rd, we added \\details explaining some implemented methods 
and outputs.

> \\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.

>Please unwrap the examples if they are executable in < 5 sec, or replace
\\dontrun{} with \\donttest{}.

Done. \\dontrun{} -> \\donttest{}

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package.
