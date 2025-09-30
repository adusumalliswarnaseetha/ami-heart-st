var spotHref = function(tabName) {
    var dropdownList = document.getElementsByTagName("a");
    for (var i = 0; i < dropdownList.length; i++) {
        var link = dropdownList[i];
        if(link.getAttribute("data-value") == tabName) {
            link.click();
        }
    }
};



var signalHref = function(tabName, geneName) {
    var dropdownList = document.getElementsByTagName("a");
    for (var i = 0; i < dropdownList.length; i++) {
        var link = dropdownList[i];
        if(link.getAttribute("data-value") == tabName) {
            link.click();
            Shiny.setInputValue('clicked_spot_gene', geneName);
        }
    }
};