function getUrlParam(p) {
    var u = window.location.search.substring(1);
    var v = u.split('&');
    for (var i = 0; i < v.length; i++) {
        var n = v[i].split('=');
        if (n[0] == p) {
	        if (n.length == 1)
	        	return true;
            return n[1];
        }
    }
    return false;
}

/**
 * If user provides "exp" GET variable then expand all showables. I.e. for development. 
 **/
if (getUrlParam("exp") !== false) {
	$(document).ready(function(){
		$('div.showable-header > a').click();
	});
}
