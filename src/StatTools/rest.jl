
using HttpServer, DataFrames


function loadScoredFile(dfd::DataFrame)
    println("loadScoredFile")
    return true
end

function rest(ip::String)
    http = HttpHandler() do req::Request, res::Response
        if ismatch(r"^/optimize/", req.resource)
            try 
                dfd = readtable(IOBuffer(convert(String,req.data) ) )
                println(dfd)
                if loadScoredFile(dfd)
                    return Response("done")
                else
                    return Response("error - could not upload file")
                end
            catch
                return Response(422, "Error: coudn't solve puzzle.")
            end    
        elseif ismatch(r"^/tst/", req.resource)
            return Response("ooo hello!")
        else
            return Response("unrecognized URL")
        end    
    end
    server = Server( http )
    run(server, host=IPv4(ip), port=8000)
end


rest("10.63.36.18")




#wget 10.63.36.18:8000/hello/name/  --post-data "email=abc@abc.com&file1=@FILE_HERE&file2=@FILE_HERE"
#wget 10.63.36.18:8000/hello/name/  --post-file=file
# wget 10.63.36.18:8000/optimize/  --post-file=/mnt/resource/analytics/models/rev/scored.csv















