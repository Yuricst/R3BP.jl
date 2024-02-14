"""
    __print_verbose(message::String, verbose::Bool)

Print message if verbose==true
"""
function __print_verbose(message::String, verbose::Bool)
    if verbose==true
        println(message)
    end
end