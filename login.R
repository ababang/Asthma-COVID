# ------------------------------------------------------------------------------
# Database access code
# ------------------------------------------------------------------------------

library(RODBC)
library(tcltk)

getlogin <- function(userName=''){
    wnd <- tktoplevel()
    user <- tclVar(userName)
    passvar <- tclVar('')

    tkgrid(tklabel(wnd,text='Username:'))
    passBox <- tkentry(wnd,textvariable = user)
    tkgrid(passBox)
    
    tkgrid(tklabel(wnd,text='Password:'))
    passBox <- tkentry(wnd,textvariable=passvar,show='*')
    tkgrid(passBox)
    
    # Hitting return will also submit password.
    tkbind(passBox, '<Return>', function() tkdestroy(wnd))
    
    # OK button.
    tkgrid(tkbutton(wnd,text='OK',command=function() tkdestroy(wnd)))
    
    # Wait for user to click OK.
    tkwait.window(wnd)
    
    password <- tclvalue(passvar)
    userName <- tclvalue(user)
    
    db <- odbcConnect('PR_SAIL', userName, password)
    return(db)
}

channel <- getlogin()



